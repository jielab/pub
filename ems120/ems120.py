import os, json, time, re, math
from pathlib import Path
from collections import Counter, defaultdict
from itertools import groupby

import pandas as pd
import torch
from transformers import (
    BertTokenizer,
    BertForSequenceClassification,
    Trainer,
    TrainingArguments,
)

# ================= PATH =================
BASE_MODEL_DIR = Path(r"D:\data\AI模型文件\hfl").as_posix()
OUTPUT_MODEL_DIR = Path(r"D:\analysis\ems120\ml\hfl").as_posix()
TRAIN_FILE = Path(r"D:\analysis\ems120\2019.train_dx.xlsx").as_posix()

TRAINED_MODEL_FILE = Path(OUTPUT_MODEL_DIR, "pytorch_model.bin").as_posix()
KEYWORDS_FILE = Path(OUTPUT_MODEL_DIR, "keywords.json").as_posix()
LABELS_FILE = Path(OUTPUT_MODEL_DIR, "labels.json").as_posix()

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
DX_TEXT_COLS = ["性别", "年龄", "呼救原因", "病因", "伤病程度", "症状", "主诉", "病史", "初步诊断", "补充诊断"]
LABEL_COL = "疾病类型"

TOKENIZER = None
MODEL = None
LABEL2ID = None
ID2LABEL = None
KEYWORD_DICT = None  # label -> {token: weight}

# ================= LOG =================
def log(msg):
    print(f"[EMS120 {time.strftime('%H:%M:%S')}] {msg}", flush=True)

def gpu_mem():
    if torch.cuda.is_available():
        log(f"GPU memory allocated: {torch.cuda.memory_allocated()/1024**3:.2f} GB")

# ================= TEXT =================
def build_texts_from_df(df: pd.DataFrame):
    cols = [c for c in DX_TEXT_COLS if c in df.columns]
    log(f"Using columns: {cols}")
    if not cols:
        return [""] * len(df)
    return df[cols].fillna("").astype(str).agg(" ".join, axis=1).tolist()

# ================= KEYWORD (TOKEN-BASED, CLEAN) =================
_re_tok = re.compile(r"[A-Za-z]+|\d+|[\u4e00-\u9fff]+")

def tokenize_text(s: str):
    """Tokenize into meaningful tokens; drop pure numbers and very short noise.
    For Chinese sequences, keep bigrams to capture phrases like '腹痛','胸痛','呼吸困难' (via bigrams).
    """
    if s is None:
        return []
    s = str(s)
    toks = []
    for m in _re_tok.finditer(s):
        t = m.group(0).strip()
        if not t:
            continue
        # drop pure digits
        if t.isdigit():
            continue
        # latin token: keep len>=2
        if re.fullmatch(r"[A-Za-z]+", t):
            if len(t) >= 2:
                toks.append(t.lower())
            continue
        # chinese token
        if re.fullmatch(r"[\u4e00-\u9fff]+", t):
            if len(t) == 1:
                continue
            if len(t) == 2:
                toks.append(t)
            else:
                # keep bigrams + full token (bounded)
                toks.append(t[:6])  # cap very long spans
                for i in range(0, min(len(t) - 1, 10)):  # avoid explosion
                    toks.append(t[i : i + 2])
            continue
    return toks

def learn_keywords(df: pd.DataFrame, topk: int = 40):
    """Learn per-label token weights; avoid punctuation/digits; return label->token->weight."""
    log("Learning keyword dictionary (token-based)...")
    texts = build_texts_from_df(df)
    labels = df[LABEL_COL].astype(str)

    global_cnt = Counter()
    label_cnt = defaultdict(Counter)

    for t, lab in zip(texts, labels):
        toks = tokenize_text(t)
        if not toks:
            continue
        c = Counter(toks)
        global_cnt.update(c)
        label_cnt[lab].update(c)

    # weight: freq_label * idf-ish (penalize very common tokens)
    # w = f_l * log(1 + N / (f_g+1))
    N = sum(global_cnt.values()) + 1
    kw = {}
    for lab, cnt in label_cnt.items():
        scores = {}
        for tok, f_l in cnt.items():
            f_g = global_cnt.get(tok, 0)
            if f_g <= 1:
                continue
            idf = math.log(1.0 + N / (f_g + 1.0))
            w = float(f_l) * idf
            # drop overly generic tokens
            if w <= 0:
                continue
            scores[tok] = w
        # keep topk
        top = sorted(scores.items(), key=lambda x: x[1], reverse=True)[:topk]
        kw[lab] = {tok: round(w, 4) for tok, w in top}

    log("Keyword learning finished")
    return kw

def keyword_predict(text: str, max_reason: int = 6):
    """Score by matched learned tokens; reason = matched top-weight tokens."""
    if not KEYWORD_DICT:
        return "其他", "NoKeywordModel"
    toks = set(tokenize_text(text))
    if not toks:
        return "其他", "EmptyText"

    best_lab = "其他"
    best_score = 0.0
    best_hits = []

    for lab, tok_w in KEYWORD_DICT.items():
        hits = [(tok, tok_w[tok]) for tok in toks if tok in tok_w]
        if not hits:
            continue
        score = sum(w for _, w in hits)
        if score > best_score:
            best_score = score
            best_lab = lab
            best_hits = sorted(hits, key=lambda x: x[1], reverse=True)

    if best_score <= 0 or not best_hits:
        return "其他", "NoMatch"

    reason = ",".join([tok for tok, _ in best_hits[:max_reason]])
    return best_lab, reason

# ================= DATASET =================
class DxDataset(torch.utils.data.Dataset):
    def __init__(self, enc, lab):
        self.enc = enc
        self.lab = lab

    def __getitem__(self, i):
        item = {k: torch.tensor(v[i]) for k, v in self.enc.items()}
        item["labels"] = torch.tensor(self.lab[i])
        return item

    def __len__(self):
        return len(self.enc["input_ids"])

def _build_label_maps_from_trainfile():
    df = pd.read_excel(TRAIN_FILE)
    if LABEL_COL not in df.columns:
        raise ValueError(f"training data missing label column: {LABEL_COL}")
    uniq = sorted(df[LABEL_COL].astype(str).unique())
    l2i = {l: i for i, l in enumerate(uniq)}
    i2l = {i: l for l, i in l2i.items()}
    return l2i, i2l

# ================= TRAIN / LOAD =================
def train_dx_model():
    global TOKENIZER, MODEL, LABEL2ID, ID2LABEL, KEYWORD_DICT
    torch.set_grad_enabled(True)

    log(f"trained model file {TRAINED_MODEL_FILE} not found")
    log(f"train model using model files {BASE_MODEL_DIR} and data file {TRAIN_FILE}")

    df = pd.read_excel(TRAIN_FILE)
    if LABEL_COL not in df.columns:
        raise ValueError(f"training data missing label column: {LABEL_COL}")

    texts = build_texts_from_df(df)
    labels = df[LABEL_COL].astype(str)

    uniq = sorted(labels.unique())
    LABEL2ID = {l: i for i, l in enumerate(uniq)}
    ID2LABEL = {i: l for l, i in LABEL2ID.items()}

    KEYWORD_DICT = learn_keywords(df)

    y = labels.map(LABEL2ID).tolist()

    log("Loading tokenizer...")
    TOKENIZER = BertTokenizer.from_pretrained(BASE_MODEL_DIR)

    log("Tokenizing training texts...")
    enc = TOKENIZER(texts, truncation=True, padding=True, max_length=128)

    log("Loading base model...")
    MODEL = BertForSequenceClassification.from_pretrained(
        BASE_MODEL_DIR, num_labels=len(uniq)
    ).to(DEVICE)

    MODEL.train()
    MODEL.requires_grad_(True)
    gpu_mem()

    args = TrainingArguments(
        output_dir=OUTPUT_MODEL_DIR,
        per_device_train_batch_size=16,
        num_train_epochs=2,
        logging_steps=50,
        logging_strategy="steps",
        disable_tqdm=False,
        save_strategy="no",
        fp16=torch.cuda.is_available(),
        report_to=[],
    )

    log("Training started...")
    Trainer(model=MODEL, args=args, train_dataset=DxDataset(enc, y)).train()

    log("Saving trained model...")
    MODEL.save_pretrained(OUTPUT_MODEL_DIR)
    TOKENIZER.save_pretrained(OUTPUT_MODEL_DIR)

    with open(KEYWORDS_FILE, "w", encoding="utf-8") as f:
        json.dump(KEYWORD_DICT, f, ensure_ascii=False)

    with open(LABELS_FILE, "w", encoding="utf-8") as f:
        json.dump(
            {"label2id": LABEL2ID, "id2label": {str(k): v for k, v in ID2LABEL.items()}},
            f,
            ensure_ascii=False,
        )

    log(f"trained model saved: {TRAINED_MODEL_FILE}")
    log("training finished")

def load_model_once():
    global TOKENIZER, MODEL, KEYWORD_DICT, LABEL2ID, ID2LABEL
    if MODEL is not None:
        return

    if not os.path.exists(TRAINED_MODEL_FILE):
        train_dx_model()
    else:
        log(f"trained model file {TRAINED_MODEL_FILE} found")

    log("Loading tokenizer/model into memory...")
    TOKENIZER = BertTokenizer.from_pretrained(OUTPUT_MODEL_DIR)

    MODEL = BertForSequenceClassification.from_pretrained(
        OUTPUT_MODEL_DIR, dtype=torch.float16
    ).to(DEVICE)

    MODEL.eval()
    MODEL.requires_grad_(False)
    torch.set_grad_enabled(False)

    # keywords
    if os.path.exists(KEYWORDS_FILE):
        log(f"keyword file {KEYWORDS_FILE} found; loading")
        with open(KEYWORDS_FILE, "r", encoding="utf-8") as f:
            KEYWORD_DICT = json.load(f)
    else:
        log(f"keyword file {KEYWORDS_FILE} not found; rebuilding from {TRAIN_FILE}")
        df = pd.read_excel(TRAIN_FILE)
        KEYWORD_DICT = learn_keywords(df)
        with open(KEYWORDS_FILE, "w", encoding="utf-8") as f:
            json.dump(KEYWORD_DICT, f, ensure_ascii=False)
        log("keyword file recreated")

    # labels
    if os.path.exists(LABELS_FILE):
        log(f"label file {LABELS_FILE} found; loading")
        with open(LABELS_FILE, "r", encoding="utf-8") as f:
            obj = json.load(f)
        LABEL2ID = obj["label2id"]
        ID2LABEL = {int(k): v for k, v in obj["id2label"].items()}
    else:
        log(f"label file {LABELS_FILE} not found; rebuilding from {TRAIN_FILE}")
        LABEL2ID, ID2LABEL = _build_label_maps_from_trainfile()
        with open(LABELS_FILE, "w", encoding="utf-8") as f:
            json.dump(
                {"label2id": LABEL2ID, "id2label": {str(k): v for k, v in ID2LABEL.items()}},
                f,
                ensure_ascii=False,
            )
        log("label file recreated")

    gpu_mem()

# ================= PHONE (PUT BEFORE DX) =================
# phone.sco is a 0-10 score: base + digit_luck + pattern_luck, then clamped to [min,max]
# phone.sco.reason is a compact "calculation string" to show how the score is formed.

# Digit-based luck (simple cultural weights)
RULE_SIMPLE = {
    "base": 4.0,
    "min": 0.0,
    "max": 10.0,
    # per-digit weights: contribution = count(digit) * weight
    "w": {"8": 1.0, "9": 1.0, "1": 0.5, "6": 0.5, "4": -4.0},
}

# Pattern-based luck (structured aesthetics) with a hard cap to avoid score saturation at 10
RULE_PATTERN = {
    "run_min": 3,     # 777 / 8888 ... (same digit run length threshold)
    "seq_min": 3,     # 123 / 987 ... (monotonic +/-1 length threshold)

    # weights for pattern contributions
    "w_run": 0.8,        # run contribution: len * w_run
    "w_seq": 0.6,        # sequence contribution: len * w_seq
    "w_repeat4": 3.5,    # 56785678 (repeat 4-digit block): fixed add
    "w_repeat2": 2.0,    # 23232323 (repeat 2-digit block): fixed add
    "w_mirror": 2.5,     # 23455432 (mirror): fixed add

    # cap the total pattern contribution
    "cap_pat": 6.0,
}

def _is_digit8(p: str) -> bool:
    return len(p) == 8 and p.isdigit()

def _repeat_block_score(p: str, block_len: int):
    """Detect repeated blocks: e.g., 56785678 (block_len=4) or 23232323 (block_len=2)."""
    if 8 % block_len != 0:
        return False, ""
    k = 8 // block_len
    blk = p[:block_len]
    if blk * k == p and len(set(blk)) > 1:
        return True, f"{blk}x{k}"
    return False, ""

def _mirror_score(p: str):
    """Detect mirror: 23455432 (a|a[::-1])."""
    a, b = p[:4], p[4:]
    if b == a[::-1] and len(set(a)) > 1:
        return True, f"{a}|{b}"
    return False, ""

def _seq_matches(p: str, seq_min: int):
    """Return monotonic +/-1 sequences as (substring, length, step)."""
    out = []
    start = 0
    cur_step = None
    cur_len = 1

    def flush(s, l, st):
        if l >= seq_min and st in (1, -1):
            out.append((p[s:s+l], l, st))

    for i in range(1, 8):
        diff = int(p[i]) - int(p[i - 1])
        if abs(diff) == 1:
            if cur_step is None:
                cur_step = diff
                cur_len = 2
                start = i - 1
            elif diff == cur_step:
                cur_len += 1
            else:
                flush(start, cur_len, cur_step)
                cur_step = diff
                cur_len = 2
                start = i - 1
        else:
            flush(start, cur_len, cur_step)
            cur_step = None
            cur_len = 1

    flush(start, cur_len, cur_step)
    return out

def _fmt_term(name: str, count: int, w: float):
    # e.g., "8(2*1.0)" or "4(1*-4.0)"
    return f"{name}({count}*{w:g})"

def eval_phone_batch(phone_list):
    """Batch scoring for 8-digit phone tail (string). Returns dict: {sco: [...], reason: [...]}."""
    scores, reasons = [], []

    w_digits = {str(k): float(v) for k, v in RULE_SIMPLE["w"].items()}
    base = float(RULE_SIMPLE["base"])
    min_s = float(RULE_SIMPLE["min"])
    max_s = float(RULE_SIMPLE["max"])

    w_run = float(RULE_PATTERN["w_run"])
    w_seq = float(RULE_PATTERN["w_seq"])
    cap_pat = float(RULE_PATTERN["cap_pat"])

    for p in list(phone_list):
        if pd.isna(p) or p is None:
            p = ""
        p = str(p).strip()
        if p.endswith(".0"):
            p = p[:-2]

        if not _is_digit8(p):
            scores.append(float("nan"))
            reasons.append("")
            continue

        terms = []
        raw = base
        terms.append(f"base({base:g})")

        # -------- digit luck --------
        digit_sum = 0.0
        for d, w in w_digits.items():
            c = p.count(d)
            if c:
                digit_sum += c * w
                terms.append(_fmt_term(d, c, w))
        raw += digit_sum

        # -------- pattern luck --------
        pat = 0.0

        # runs (e.g., 777 / 8888)
        for k, g in groupby(p):
            l = sum(1 for _ in g)
            if l >= int(RULE_PATTERN["run_min"]):
                add = l * w_run
                pat += add
                terms.append(f"run {k*l}({l}*{w_run:g})")

        # sequences (e.g., 123 / 987)
        for sub, l, st in _seq_matches(p, int(RULE_PATTERN["seq_min"])):
            add = l * w_seq
            pat += add
            terms.append(f"seq {sub}({l}*{w_seq:g})")

        # repeat blocks
        ok, tag = _repeat_block_score(p, 4)
        if ok:
            add = float(RULE_PATTERN["w_repeat4"])
            pat += add
            terms.append(f"rep4 {tag}({add:g})")
        ok, tag = _repeat_block_score(p, 2)
        if ok:
            add = float(RULE_PATTERN["w_repeat2"])
            pat += add
            terms.append(f"rep2 {tag}({add:g})")

        # mirror
        ok, tag = _mirror_score(p)
        if ok:
            add = float(RULE_PATTERN["w_mirror"])
            pat += add
            terms.append(f"mirror {tag}({add:g})")

        # cap pattern contribution (avoid over-saturation at 10)
        if pat > cap_pat:
            terms.append(f"pat_cap({cap_pat:g})")
            pat = cap_pat

        raw += pat

        final = max(min_s, min(max_s, round(raw, 1)))

        # compact calculation string
        reason = " + ".join(terms) + f" = {final:g}"
        scores.append(final)
        reasons.append(reason)

    return {"sco": scores, "reason": reasons}

# ================= DX =================
def eval_dx_batch(df, data_name=None):
    load_model_once()
    if data_name is None:
        data_name = f"rows={len(df)}"
    log(f"now process data {data_name}")

    texts = build_texts_from_df(df)

    kw, kw_reason = [], []
    for t in texts:
        k, r = keyword_predict(t)
        kw.append(k)
        kw_reason.append(r)

    BATCH_SIZE = 128
    preds_all, conf_all = [], []
    log(f"Starting inference: {len(texts)} samples (batch={BATCH_SIZE})")

    for i in range(0, len(texts), BATCH_SIZE):
        if i % 2000 == 0:
            log(f"Inference progress: {i}/{len(texts)}")

        batch = texts[i : i + BATCH_SIZE]
        enc = TOKENIZER(batch, truncation=True, padding=True, max_length=128, return_tensors="pt")
        enc = {k: v.to(DEVICE) for k, v in enc.items()}

        with torch.no_grad(), torch.amp.autocast("cuda", dtype=torch.float16):
            logits = MODEL(**enc).logits
            probs = torch.softmax(logits, dim=1)

        preds = torch.argmax(probs, dim=1).cpu().numpy()
        conf = probs.max(dim=1).values.cpu().numpy()

        preds_all.extend(preds)
        conf_all.extend(conf)

        del enc, logits, probs
        torch.cuda.empty_cache()

    ml = [ID2LABEL[int(i)] for i in preds_all]
    ml_reason = [f"{ID2LABEL[int(i)]} (p={c:.2f})" for i, c in zip(preds_all, conf_all)]

    log("Inference finished")
    return {"kw": kw, "kw_reason": kw_reason, "ml": ml, "ml_reason": ml_reason}


    

# ================= 房价地图 =================    
# import pandas as pd
# import keplergl
# dat = pd.read_csv("/mnt/d/files/120.txt", sep = "\t"); # dat.head()
# map = keplergl.KeplerGl(height = 500); map # 需要把这个作为最后一行命令
# map.add_data(data = dat.copy(), name = "house")
# %run 120.config.py
# with open('120.config.py', 'w') as f: f.write('config = {}'.format(map.config))