"""Launch a supervised fine-tuning job on OpenAI or Azure OpenAI.

Defaults reproduce the manuscript configuration (Key resources table): base model
gpt-4.1-2025-04-14, batch size 8, learning-rate multiplier 1.2, 3 epochs, trained on
``data/training/training_v3.jsonl.gz``. The resulting model is private to the account that
trains it; pass its deployment name to ``pathwayseeker ask --model``.

    pathwayseeker finetune data/training/training_v3.jsonl.gz --provider azure
"""

import gzip
import shutil
import tempfile
from pathlib import Path
from typing import Optional

PAPER_CONFIG = {"model": "gpt-4.1-2025-04-14", "n_epochs": 3, "batch_size": 8,
                "learning_rate_multiplier": 1.2}


def _client(provider: str):
    if provider == "azure":
        import os

        from openai import AzureOpenAI

        return AzureOpenAI(azure_endpoint=os.environ.get("AZURE_OPENAI_ENDPOINT"),
                           api_key=os.environ.get("AZURE_OPENAI_API_KEY"),
                           api_version=os.environ.get("AZURE_OPENAI_API_VERSION", "2025-01-01-preview"))
    from openai import OpenAI

    return OpenAI()


def _plain_jsonl(path: Path) -> Path:
    if path.suffix != ".gz":
        return path
    tmp = Path(tempfile.mkdtemp()) / path.stem
    with gzip.open(path, "rb") as src, open(tmp, "wb") as dst:
        shutil.copyfileobj(src, dst)
    return tmp


def create_job(training_file: str, provider: str = "azure", model: str = PAPER_CONFIG["model"],
               n_epochs: int = PAPER_CONFIG["n_epochs"], batch_size: int = PAPER_CONFIG["batch_size"],
               learning_rate_multiplier: float = PAPER_CONFIG["learning_rate_multiplier"],
               validation_file: Optional[str] = None, suffix: str = "pathwayseeker"):
    """Upload the training (and optional validation) JSONL and create the fine-tuning job."""
    client = _client(provider)

    def upload(p):
        with open(_plain_jsonl(Path(p)), "rb") as fh:
            return client.files.create(file=fh, purpose="fine-tune").id

    kwargs = {"validation_file": upload(validation_file)} if validation_file else {}
    return client.fine_tuning.jobs.create(
        training_file=upload(training_file), model=model, suffix=suffix,
        hyperparameters={"n_epochs": n_epochs, "batch_size": batch_size,
                         "learning_rate_multiplier": learning_rate_multiplier},
        **kwargs)
