from .dvh import DVH


def homogeneity_index(dvh: DVH, target_dose: float) -> float:
    return (dvh.D(0.05) - dvh.D(0.95)) / target_dose
