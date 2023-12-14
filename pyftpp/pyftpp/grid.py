from dataclasses import dataclass
from typing import Tuple
import numpy as np


@dataclass
class BoundingBox:
    x_min: float
    x_max: float
    y_min: float
    y_max: float
    z_min: float
    z_max: float


@dataclass
class Grid:
    shape: Tuple[int, int, int]
    # Spacing in mm. (?)
    spacing: Tuple[float, float, float]
    # Origin coordinates in mm. (?)
    origin: Tuple[float, float, float]

    @property
    def bounding_box(self) -> BoundingBox:
        assert (self.spacing[0] >= 0) \
            and (self.spacing[1] >= 0) \
            and (self.spacing[2] >= 0)

        return BoundingBox(
            self.origin[0],
            self.origin[0] + self.shape[0]*self.spacing[0],
            self.origin[1],
            self.origin[1] + self.shape[1]*self.spacing[1],
            self.origin[2],
            self.origin[2] + self.shape[2]*self.spacing[2],
        )

    def is_finer_than(self, other) -> bool:
        assert np.isclose(self.spacing[0], self.spacing[1], atol=1e-6)
        assert np.isclose(other.spacing[0], other.spacing[1], atol=1e-6)
        assert np.isclose(self.spacing[2], other.spacing[2], atol=1e-6)

        return self.spacing[0] <= other.spacing[0]

    def __eq__(self, other) -> bool:
        atol = 1e-6

        if not isinstance(other, Grid):
            raise NotImplemented

        for i in range(3):
            if self.shape[i] != other.shape[i]:
                return False
            if abs(self.origin[i] - other.origin[i]) > atol:
                return False
            if abs(self.spacing[i] - other.spacing[i]) > atol:
                return False
        return True
