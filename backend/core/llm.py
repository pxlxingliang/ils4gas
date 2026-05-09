from abc import ABC, abstractmethod
from typing import Iterator, List, Dict, Any, Optional, AsyncIterator

from openai import OpenAI, AsyncOpenAI


class LLMProvider(ABC):
    @abstractmethod
    def invoke(self, messages: List[Dict], **kwargs) -> str:
        ...

    @abstractmethod
    def stream_invoke(self, messages: List[Dict], **kwargs) -> Iterator[str]:
        ...

    @abstractmethod
    async def ainvoke(self, messages: List[Dict], **kwargs) -> str:
        ...

    @abstractmethod
    async def astream_invoke(
        self, messages: List[Dict], **kwargs
    ) -> AsyncIterator[str]:
        ...

    @property
    @abstractmethod
    def model_name(self) -> str:
        ...

    @property
    @abstractmethod
    def context_limit(self) -> int:
        ...


class OpenAICompatibleProvider(LLMProvider):
    def __init__(
        self,
        api_key: str,
        base_url: Optional[str],
        model_name: str,
        context_limit: int = 128000,
    ):
        self._model_name = model_name
        self._context_limit = context_limit
        self.client = OpenAI(api_key=api_key, base_url=base_url)
        self.async_client = AsyncOpenAI(api_key=api_key, base_url=base_url)

    def invoke(self, messages: List[Dict], **kwargs) -> str:
        resp = self.client.chat.completions.create(
            model=self._model_name, messages=messages, **kwargs
        )
        return resp.choices[0].message.content or ""

    def stream_invoke(self, messages: List[Dict], **kwargs) -> Iterator[str]:
        stream = self.client.chat.completions.create(
            model=self._model_name, messages=messages, stream=True, **kwargs
        )
        for chunk in stream:
            delta = chunk.choices[0].delta
            if delta.content:
                yield delta.content

    async def ainvoke(self, messages: List[Dict], **kwargs) -> str:
        resp = await self.async_client.chat.completions.create(
            model=self._model_name, messages=messages, **kwargs
        )
        return resp.choices[0].message.content or ""

    async def astream_invoke(
        self, messages: List[Dict], **kwargs
    ) -> AsyncIterator[str]:
        stream = await self.async_client.chat.completions.create(
            model=self._model_name, messages=messages, stream=True, **kwargs
        )
        async for chunk in stream:
            delta = chunk.choices[0].delta
            if delta.content:
                yield delta.content

    @property
    def model_name(self) -> str:
        return self._model_name

    @property
    def context_limit(self) -> int:
        return self._context_limit
