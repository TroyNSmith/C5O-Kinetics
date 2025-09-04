from IPython.display import Javascript, display
from ipywidgets import (
    Button,
    Combobox,
    Dropdown,
    IntText,
    Label,
    Layout,
    RadioButtons,
    Tab,
    Text,
)


# --- Utility: Copy to clipboard ---
def copy_to_clipboard(value):
    display(Javascript(f"navigator.clipboard.writeText('{value}')"))


# --- Widget Factory Functions ---
def combobox(options, description):
    return Combobox(
        options=options,
        description=description,
        layout=Layout(height="25", width="auto"),
        ensure_option=False,
    )


def dropdown(options, description=""):
    if type(options[0]) is tuple:
        default = options[0][1]
    else:
        default = options[0]
    return Dropdown(
        options=options,
        value=default,
        description=description,
        layout=Layout(height="25", width="auto"),
    )


def button(description, style="info", **kwargs):
    button = Button(
        description=description,
        layout=Layout(height="auto", width="auto"),
        **kwargs,
    )
    button.style.button_color = "#AAAAAA"
    return button


def label(value):
    label = Label(value=value, layout=Layout(height="auto", width="auto"))
    return label


def radio(description, options):
    radio_buttons = RadioButtons(
        options=options,
        layout={"width": "max-content", "height": "100%"},
        description=description,
        disabled=False,
    )
    return radio_buttons


def text(placeholder: str, value: str, description: str = ""):
    return Text(
        placeholder=placeholder,
        value=value,
        description=description,
        layout=Layout(height="30%", width="auto"),
    )


def int_text(value: int = 2, step: int = 1, description: str = ""):
    return IntText(
        value=value,
        description=description,
    )


def tab(children: list, titles: list):
    return Tab(children=children, titles=titles)
