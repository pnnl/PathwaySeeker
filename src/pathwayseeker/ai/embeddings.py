"""Embedding store using Azure OpenAI (lazy initialization)."""

import os
import pickle

import numpy as np
from tqdm.auto import tqdm

from . import config as C


def _get_client():
    """Lazily create Azure OpenAI client (only when actually needed)."""
    from openai import AzureOpenAI

    endpoint = os.getenv("ENDPOINT_URL", "https://pathways.openai.azure.com/")
    api_key = os.getenv("AZURE_OPENAI_API_KEY_OMICS")
    deployment = os.getenv("DEPLOYMENT_NAME", "text-embedding-3-large")
    api_version = "2024-12-01-preview"

    return AzureOpenAI(
        azure_endpoint=endpoint,
        api_key=api_key,
        api_version=api_version,
    )


def _batched(iterable, n):
    """Yield successive n-sized chunks from iterable."""
    it = iter(iterable)
    while True:
        chunk = []
        try:
            for _ in range(n):
                chunk.append(next(it))
        except StopIteration:
            if chunk:
                yield chunk
            break
        if chunk:
            yield chunk


class EmbeddingStore:
    """
    Lazily builds / loads OpenAI embeddings for node textual
    representations (enzyme description, compound name, etc.).
    """

    def __init__(self, node_text_lookup: dict[str, str]):
        self.lookup = node_text_lookup
        self.vectors = self._load_or_build()

    def _load_or_build(self) -> dict[str, np.ndarray]:
        if os.path.exists(C.EMBED_CACHE_FILE):
            with open(C.EMBED_CACHE_FILE, "rb") as fh:
                return pickle.load(fh)

        client = _get_client()
        vecs: dict[str, np.ndarray] = {}
        for chunk in _batched(self.lookup.items(), C.EMBED_BATCH):
            ids, texts = zip(*chunk)
            resp = client.embeddings.create(
                model=C.EMBED_MODEL,
                input=list(texts),
            )
            for _id, emb in zip(ids, resp.data):
                vecs[_id] = np.asarray(emb.embedding, dtype=np.float32)

        os.makedirs(os.path.dirname(C.EMBED_CACHE_FILE), exist_ok=True)
        with open(C.EMBED_CACHE_FILE, "wb") as fh:
            pickle.dump(vecs, fh)
        return vecs

    def vector(self, node_id: str) -> np.ndarray:
        return self.vectors[node_id]

    def similarity(self, id1: str, id2: str, metric: str = C.LP_SIM_METRIC) -> float:
        v1, v2 = self.vector(id1), self.vector(id2)
        if metric == "cosine":
            return float(v1 @ v2 / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-9))
        return float(v1 @ v2)
