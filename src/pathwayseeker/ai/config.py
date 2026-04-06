"""AI configuration constants."""

# Model & embedding settings
EMBED_MODEL = "text-embedding-3-small"
EMBED_BATCH = 96
EMBED_CACHE_FILE = "cache/embeddings.pkl"

# Link-prediction settings
LP_SIM_METRIC = "cosine"
LP_TOP_K = 10
LP_MIN_SIM = 0.32
LP_WEIGHT_SCALE = 10.0

# Search variants
VARIANTS = {
    "baseline":            {"use_embed": False, "use_llm": False, "use_lp": False},
    "use_embedding":       {"use_embed": True,  "use_llm": False, "use_lp": False},
    "use_llm":             {"use_embed": True,  "use_llm": True,  "use_lp": False},
    "use_link_prediction": {"use_embed": True,  "use_llm": True,  "use_lp": True},
}
