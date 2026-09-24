"""Minimal provider-agnostic JSON chat interface for the reasoning layer.

Any object with a ``complete_json(prompt, system=None, temperature=0.3) -> dict`` method works.
Built-in adapters cover OpenAI (and OpenAI-compatible servers such as vLLM or Ollama through
``base_url``), Azure OpenAI, and Anthropic.

Environment variables used by :func:`get_llm`:

``PATHWAYSEEKER_LLM``      provider: ``openai`` (default), ``azure`` or ``anthropic``
``PATHWAYSEEKER_MODEL``    model or Azure deployment name
``OPENAI_API_KEY``, ``OPENAI_BASE_URL``
``AZURE_OPENAI_API_KEY``, ``AZURE_OPENAI_ENDPOINT``, ``AZURE_OPENAI_API_VERSION``
``ANTHROPIC_API_KEY``
"""

from __future__ import annotations

import json
import os
import re
from typing import Optional, Protocol

DEFAULT_MODELS = {"openai": "gpt-4.1", "azure": "gpt-4.1", "anthropic": "claude-sonnet-5"}


class LLM(Protocol):
    def complete_json(self, prompt: str, system: Optional[str] = None,
                      temperature: float = 0.3) -> dict: ...


def parse_json(text: str) -> dict:
    """Parse the first JSON object in ``text``; return {} if none parses."""
    if not text:
        return {}
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        pass
    fenced = re.search(r"```(?:json)?\s*(\{.*?\})\s*```", text, re.S)
    if fenced:
        try:
            return json.loads(fenced.group(1))
        except json.JSONDecodeError:
            pass
    start = text.find("{")
    end = text.rfind("}")
    if start != -1 and end > start:
        try:
            return json.loads(text[start:end + 1])
        except json.JSONDecodeError:
            return {}
    return {}


class OpenAILLM:
    """OpenAI chat completions in JSON mode. Pass ``base_url`` for OpenAI-compatible servers."""

    def __init__(self, model: str = DEFAULT_MODELS["openai"], client=None, **client_kw):
        if client is None:
            from openai import OpenAI

            client = OpenAI(**client_kw)
        self.client = client
        self.model = model

    def complete_json(self, prompt, system=None, temperature=0.3):
        messages = ([{"role": "system", "content": system}] if system else []) + [
            {"role": "user", "content": prompt}]
        resp = self.client.chat.completions.create(
            model=self.model, messages=messages, temperature=temperature,
            response_format={"type": "json_object"})
        return parse_json(resp.choices[0].message.content)


class AzureOpenAILLM(OpenAILLM):
    """Azure OpenAI; ``model`` is the deployment name (fine-tuned deployments included)."""

    def __init__(self, model: str = DEFAULT_MODELS["azure"], **client_kw):
        from openai import AzureOpenAI

        client_kw.setdefault("azure_endpoint", os.environ.get("AZURE_OPENAI_ENDPOINT"))
        client_kw.setdefault("api_key", os.environ.get("AZURE_OPENAI_API_KEY"))
        client_kw.setdefault("api_version", os.environ.get("AZURE_OPENAI_API_VERSION", "2024-12-01-preview"))
        super().__init__(model=model, client=AzureOpenAI(**client_kw))


class AnthropicLLM:
    """Anthropic Messages API; the JSON object is parsed from the text reply."""

    def __init__(self, model: str = DEFAULT_MODELS["anthropic"], client=None, max_tokens: int = 4096,
                 **client_kw):
        if client is None:
            import anthropic

            client = anthropic.Anthropic(**client_kw)
        self.client = client
        self.model = model
        self.max_tokens = max_tokens

    def complete_json(self, prompt, system=None, temperature=0.3):
        kw = {"system": system} if system else {}
        resp = self.client.messages.create(
            model=self.model, max_tokens=self.max_tokens, temperature=temperature,
            messages=[{"role": "user", "content": prompt + "\n\nRespond with one JSON object and nothing else."}],
            **kw)
        text = "".join(b.text for b in resp.content if getattr(b, "type", "") == "text")
        return parse_json(text)


def get_llm(provider: Optional[str] = None, model: Optional[str] = None) -> LLM:
    """Build an LLM adapter from arguments or environment variables."""
    provider = (provider or os.environ.get("PATHWAYSEEKER_LLM", "openai")).lower()
    model = model or os.environ.get("PATHWAYSEEKER_MODEL") or DEFAULT_MODELS.get(provider)
    if provider == "openai":
        return OpenAILLM(model=model)
    if provider == "azure":
        return AzureOpenAILLM(model=model)
    if provider == "anthropic":
        return AnthropicLLM(model=model)
    raise ValueError(f"Unknown LLM provider {provider!r}; use openai, azure or anthropic")
