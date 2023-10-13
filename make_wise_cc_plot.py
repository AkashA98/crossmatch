import copy, matplotlib
from svgpath2mpl import parse_path
from dataclasses import dataclass
from typing import Dict, Any, Tuple
from matplotlib import pyplot as plt
from matplotlib.patches import PathPatch
from matplotlib.transforms import Affine2D


default_wise_objects = [
    "COOL_T_DWARFS",
    "STARS",
    "ELLIPTICALS",
    "SPIRALS",
    "LIRGS",
    "STARBURSTS",
    "SEYFERTS",
    "QSOS",
    "OBSCURED_AGN",
]


class WiseClass:
    """WISE object classes defined in the WISE color-color plot."""

    COOL_T_DWARFS = "CoolTDwarfs"
    STARS = "Stars"
    ELLIPTICALS = "Ellipticals"
    SPIRALS = "Spirals"
    LIRGS = "LIRGs"
    STARBURSTS = "Starbursts"
    SEYFERTS = "Seyferts"
    QSOS = "QSOs"
    OBSCURED_AGN = "ObscuredAGN"


class Wisesvgpaths:
    """Class object that stores svg paths for WISE objects"""

    WISE_MPL_PATHS = {
        WiseClass.COOL_T_DWARFS: "M2.19,2.29C2.64,2.29,3,1.29,3,0H1.3C1.34,1.3,1.74,2.29,2.19,2.29Z",
        WiseClass.LIRGS: "M3.58,3.22a.91.91,0,0,1,1.06-1A2.82,2.82,0,0,1,6.55,3.31c.38.56-.1,1-1.06,1S3.67,3.8,3.58,3.22",
        WiseClass.OBSCURED_AGN: "M6.36,2.32a4.78,4.78,0,0,0-1.43-.8,2.17,2.17,0,0,0-.5-.11,1.78,1.78,0,0,0-.34,0c-.37,0-.62.06-.64.14s.3.26.72.38L4.35,2a3.2,3.2,0,0,1,.78.52c.62.59,1,1.07,1.37,1S7,2.94,6.36,2.32Z",
        WiseClass.STARS: "M1.67,3.39H1.58c-.36,0-.64.17-.64.4s.21.33.51.33l.2,0c.38-.08.69-.25.7-.39S2.06,3.42,1.67,3.39Z",
        WiseClass.ELLIPTICALS: "M2.33,3.82A.5.5,0,0,0,2,3.69a.64.64,0,0,0-.2,0c-.22,0-.39.14-.39.2s.26.12.59.13H2A.6.6,0,0,0,2.33,4,.09.09,0,0,0,2.33,3.82Z",
        WiseClass.STARBURSTS: "M5.34,2.74A.68.68,0,0,0,5,2.65c-.31,0-.5.26-.5.67s.22.74.56.74a.69.69,0,0,0,.24,0c.42-.17.76-.43.77-.58S5.76,3,5.34,2.74Z",
        WiseClass.SPIRALS: "M5.29,3.5c-.09-.11-.37-.17-.8-.17a7.1,7.1,0,0,0-1.07.1A4.07,4.07,0,0,0,2,3.78a.06.06,0,0,0,0,.09,3.39,3.39,0,0,0,1.51.22l.47,0A2.55,2.55,0,0,0,5.26,3.7C5.33,3.63,5.34,3.57,5.29,3.5Z",
        WiseClass.SEYFERTS: "M3.28,2.86c0-.32.48-.57,1.11-.57s1.15.26,1.15.59-.5.57-1.11.56-1.12-.27-1.15-.58",
        WiseClass.QSOS: "M4.56,2.62a2,2,0,0,1-.81.66c-.22.05-.36-.09-.3-.33a1.22,1.22,0,0,1,.81-.66c.4-.13.53,0,.3.33",
    }


@dataclass
class WisePatchConfig:
    """Style and annotation configurations for the patch drawn to represent a
    WISE object class in the WISE color-color plot.

    Attributes:
        style (Dict[str, Any]): Any style keyword arguments and values
            supported by `matplotlib.patches.PathPatch`.
        annotation_text (str): Text to annotate the patch.
        annotation_position (Tuple[float, float]): Position in data coordinates
        for the annotation text.
    """

    style: Dict[str, Any]
    annotation_text: str
    annotation_position: Tuple[float, float]

    def copy(self):
        return copy.deepcopy(self)


WISE_DEFAULT_PATCH_CONFIGS = {
    WiseClass.COOL_T_DWARFS: WisePatchConfig(
        style=dict(fc="#cb4627", ec="none"),
        annotation_text="Cool\nT-Dwarfs",
        annotation_position=(1.15, 3.0),
    ),
    WiseClass.STARS: WisePatchConfig(
        style=dict(fc="#e8e615", ec="none"),
        annotation_text="Stars",
        annotation_position=(0.5, 0.4),
    ),
    WiseClass.ELLIPTICALS: WisePatchConfig(
        style=dict(fc="#95c53d", ec="none"),
        annotation_text="Ellipticals",
        annotation_position=(1.0, -0.25),
    ),
    WiseClass.SPIRALS: WisePatchConfig(
        style=dict(fc="#bbdeb5", ec="none", alpha=0.7),
        annotation_text="Spirals",
        annotation_position=(2.5, 0.35),
    ),
    WiseClass.LIRGS: WisePatchConfig(
        style=dict(fc="#ecc384", ec="none"),
        annotation_text="LIRGs",
        annotation_position=(5.0, -0.1),
    ),
    WiseClass.STARBURSTS: WisePatchConfig(
        style=dict(fc="#e8e615", ec="none", alpha=0.7),
        annotation_text="ULIRGs\nLINERs\nStarbursts",
        annotation_position=(4.75, 0.5),
    ),
    WiseClass.SEYFERTS: WisePatchConfig(
        style=dict(fc="#45c7f0", ec="none", alpha=0.7),
        annotation_text="Seyferts",
        annotation_position=(3.5, 0.9),
    ),
    WiseClass.QSOS: WisePatchConfig(
        style=dict(fc="#b4e2ec", ec="none"),
        annotation_text="QSOs",
        annotation_position=(3.1, 1.25),
    ),
    WiseClass.OBSCURED_AGN: WisePatchConfig(
        style=dict(fc="#faa719", ec="none", alpha=0.8),
        annotation_text="ULIRGs/LINERs\nObscured AGN",
        annotation_position=(4.5, 1.75),
    ),
}


def get_mpl_paths(wise_obj):
    """Helper function to transform the svg paths to matplotlib paths

    Args:
        wise_obj (str): The name of the WISE object
    """
    # Define a transform
    mpl_transform = Affine2D().scale(sx=1, sy=-1).translate(-1, 4)
    p = Wisesvgpaths.WISE_MPL_PATHS[wise_obj]
    mpl_path = parse_path(p).transformed(transform=mpl_transform)

    return mpl_path


def wise_color_color_plot(wise_objects=default_wise_objects):
    """Make an empty WISE color-color plot with common object classes drawn as
    patches. The patches have default styles that may be overridden.

    Returns:
        `matplotlib.figure.Figure`: the WISE color-color figure. Access the
            axes with the `.axes` attribute.
    """
    # set the WISE object classification patch styles
    patch_styles = WISE_DEFAULT_PATCH_CONFIGS

    fig, ax = plt.subplots(figsize=(8, 6))
    for w in wise_objects:
        name = getattr(WiseClass, w)
        patch_style = patch_styles[name]
        patch_mpl = get_mpl_paths(name)
        patch = PathPatch(patch_mpl, **patch_style.style)
        ax.add_patch(patch)
        ax.annotate(
            patch_style.annotation_text,
            patch_style.annotation_position,
            ha="center",
            fontsize="medium",
        )
    ax.set_xlim(-1, 6)
    ax.set_ylim(-0.5, 4)
    ax.set_aspect(1)
    ax.set_xlabel("[4.6] - [12] (mag)", fontsize=15)
    ax.set_ylabel("[3.4] - [4.6] (mag)", fontsize=15)
    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(1))
    ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(1))
    ax.tick_params(labelsize=15)
    plt.tight_layout()
    return fig, ax


if __name__ == "__main__":
    fig, ax = wise_color_color_plot(default_wise_objects)
    plt.show()
