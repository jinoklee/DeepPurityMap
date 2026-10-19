# =============================================================================
# 4_predict.py   Track 1 step 4 purity (one sample)
#
#   3_features/{id}.seg.png → purity model → 4_purity/{id}.purity.tsv
#
#   Usage
#     python 4_predict.py --sample ID --features-dir RESULTS/ID/3_features --out RESULTS/ID/4_purity
#     python 4_predict.py --sample ID --image path/to/ID.seg.png --out ...
#
#   The model is described by a model card (default: resources/model/model.json, independent of the
#   working directory): name · version · weights file + sha256 · backbone + its structure config.
#   The weights' sha256 and backbone are verified before predicting; name · version · sha256 go into the output.
#   ImageNet pretrained weights are not needed — the checkpoint overwrites everything, so only the structure is built.
#
#   --device auto (default): cuda:0 if a GPU is available, otherwise CPU.
# =============================================================================

import argparse
import hashlib
import json
import sys
from pathlib import Path

import torch
import torch.nn as nn
from PIL import Image
from torchvision import transforms

ROOT = Path(__file__).resolve().parent.parent            # repository root
sys.path.insert(0, str(ROOT / "lib"))
from dinat import DinatConfig, DinatModel           # noqa: E402  (lib/dinat)

DEFAULT_CARD = ROOT / "resources" / "model" / "model.json"

# Used only when the checkpoint has no preprocess entry (training-time values)
FALLBACK_PREPROCESS = {
    "input_size": (224, 224),
    "norm_mean": (0.03748157247900963, 0.06450977921485901, 0.07612253725528717),
    "norm_std":  (0.09910687059164047, 0.1578933149576187,  0.17405302822589874),
}


# ── Model (same structure · names as the training code — state_dict keys are backbone.* / head.*) ──
def get_simple_classifier(feature_dim):
    hidden_dim = 128 if feature_dim > 1000 else 64 if feature_dim >= 512 else 32
    return nn.Sequential(
        nn.Dropout(0.5),
        nn.Linear(feature_dim, hidden_dim),
        nn.ReLU(),
        nn.Dropout(0.3),
        nn.Linear(hidden_dim, 1),
        nn.Sigmoid(),
    )


class DiNATRegressor(nn.Module):
    def __init__(self, config: DinatConfig):
        super().__init__()
        self.backbone = DinatModel(config)
        feature_dim = int(config.embed_dim * 2 ** (len(config.depths) - 1))   # base: 1024
        self.head = get_simple_classifier(feature_dim)

    def forward(self, x):
        out = self.backbone(x)
        feat = out.last_hidden_state.mean(dim=[1, 2])          # (B,7,7,C) → (B,C)
        return self.head(feat)


# ── Model card ───────────────────────────────────────────────────────────────
def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_card(card_path: Path) -> dict:
    if not card_path.is_file():
        raise SystemExit(f"model card not found: {card_path}")
    card = json.loads(card_path.read_text(encoding="utf-8"))
    for key in ("name", "version", "backbone", "weights", "backbone_config"):
        if not card.get(key):
            raise SystemExit(f"model card {card_path.name} is missing '{key}'")
    base = card_path.parent                                  # files are relative to the card
    card["weights_path"] = base / card["weights"]
    card["config_path"] = base / card["backbone_config"]
    for k, what in (("weights_path", "model weights"), ("config_path", "backbone structure config")):
        if not card[k].is_file():
            raise SystemExit(f"{what} not found: {card[k]}")
    actual = sha256(card["weights_path"])
    if card.get("weights_sha256") and actual != card["weights_sha256"]:
        raise SystemExit(f"weights sha256 mismatch for {card['weights']} — expected {card['weights_sha256'][:12]}…, "
                         f"got {actual[:12]}… → the weights file is not the one this model card describes")
    card["weights_sha256"] = actual
    return card


def pick_device(want: str) -> torch.device:
    if want == "cpu":
        return torch.device("cpu")
    if not torch.cuda.is_available():
        if want != "auto":
            print(f"warning: GPU not available → running on CPU (requested: {want})", file=sys.stderr)
        return torch.device("cpu")
    if want in ("auto", "cuda"):
        return torch.device("cuda:0")
    idx = int(want.split(":")[1])
    if idx >= torch.cuda.device_count():
        print(f"warning: {want} not found ({torch.cuda.device_count()} GPUs) → cuda:0", file=sys.stderr)
        return torch.device("cuda:0")
    return torch.device(want)


def load_model(card: dict, device: torch.device):
    ckpt = torch.load(card["weights_path"], map_location=device, weights_only=True)
    arch = ckpt.get("model_arch", {}).get("backbone_name", card["backbone"])
    if arch != card["backbone"]:
        raise SystemExit(f"checkpoint backbone is {arch}, but the model card says {card['backbone']}")
    model = DiNATRegressor(DinatConfig.from_json_file(str(card["config_path"])))
    model.load_state_dict(ckpt["model_state_dict"])            # strict: fails if any key differs
    model.to(device).eval()
    return model, ckpt.get("preprocess") or FALLBACK_PREPROCESS


def predict(model, preprocess, img_path: Path, device) -> float:
    tf = transforms.Compose([
        transforms.Resize(tuple(preprocess["input_size"])),
        transforms.ToTensor(),
        transforms.Normalize(mean=list(preprocess["norm_mean"]), std=list(preprocess["norm_std"])),
    ])
    x = tf(Image.open(img_path).convert("RGB")).unsqueeze(0).to(device)
    with torch.no_grad():
        return model(x).squeeze().item()


def main(argv=None):
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, ValueError):
            pass
    ap = argparse.ArgumentParser(description="Tumor purity prediction from a seg image (one sample)")
    ap.add_argument("--sample", required=True, help="sample ID")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--features-dir", help="3_features folder (looks for {sample}.seg.png)")
    src.add_argument("--image", help="path to seg.png directly")
    ap.add_argument("--out", required=True, help="output folder (4_purity)")
    ap.add_argument("--model-card", default=str(DEFAULT_CARD), help="model card JSON (default: resources/model/model.json)")
    ap.add_argument("--device", default="auto", help="auto | cpu | cuda | cuda:N")
    a = ap.parse_args(argv)

    img = Path(a.image) if a.image else Path(a.features_dir) / f"{a.sample}.seg.png"
    if not img.is_file():
        raise SystemExit(f"seg image not found: {img}")

    card = load_card(Path(a.model_card))
    device = pick_device(a.device)
    model, preprocess = load_model(card, device)
    purity = predict(model, preprocess, img, device)

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    out_file = out / f"{a.sample}.purity.tsv"
    cols = ["sample_id", "predicted_purity", "model", "model_version", "device", "weights_sha256"]
    vals = [a.sample, f"{purity:.4f}", card["name"], card["version"], device.type, card["weights_sha256"][:12]]
    with open(out_file, "w", encoding="utf-8") as f:
        f.write("\t".join(cols) + "\n" + "\t".join(vals) + "\n")
    print(f"{a.sample}  purity={purity:.4f}  ({card['name']} v{card['version']}, {device})  → {out_file}")


if __name__ == "__main__":
    main()
