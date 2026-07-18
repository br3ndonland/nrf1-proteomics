from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Protocol, cast

from adjustText import adjust_text as _untyped_adjust_text
from matplotlib.artist import Artist
from matplotlib.axes import Axes
from matplotlib.text import Text


class _AdjustText(Protocol):
    def __call__(
        self,
        texts: Sequence[Text],
        *,
        x: Sequence[float],
        y: Sequence[float],
        objects: Sequence[Artist] | None,
        ax: Axes,
        iter_lim: int,
        min_arrow_len: float,
        arrowprops: Mapping[str, object],
    ) -> object: ...


_adjust_text = cast(_AdjustText, _untyped_adjust_text)


def adjust_text_labels(
    labels: Sequence[Text],
    *,
    x: Sequence[float],
    y: Sequence[float],
    objects: Sequence[Artist] | None,
    ax: Axes,
    iter_lim: int,
    min_arrow_len: float,
    arrowprops: Mapping[str, object],
) -> None:
    """Adjust label positions through the typed subset of adjustText we use."""

    _ = _adjust_text(
        labels,
        x=x,
        y=y,
        objects=objects,
        ax=ax,
        iter_lim=iter_lim,
        min_arrow_len=min_arrow_len,
        arrowprops=arrowprops,
    )
