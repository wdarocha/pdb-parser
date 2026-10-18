"""Shared type aliases for remapping layout, ordering, and clique helpers."""

from __future__ import annotations

from typing import TypeAlias

IntervalBounds: TypeAlias = tuple[float, float]
IntervalLookup: TypeAlias = dict[tuple[int, int], IntervalBounds]
AtomRef: TypeAlias = tuple[int, str]
AtomRefList: TypeAlias = list[AtomRef]
AtomLookup: TypeAlias = dict[tuple[int, str], int]
CliqueAtomSlot: TypeAlias = AtomRef | None
CliqueAtomRow: TypeAlias = list[CliqueAtomSlot]
CliqueMatrixEntry: TypeAlias = tuple[CliqueAtomRow, int]
RemapAtomLayout: TypeAlias = str
