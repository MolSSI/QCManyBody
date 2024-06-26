import re
import inspect
import importlib
from enum import Enum
from collections import namedtuple
from typing import List, Dict, Optional

import tabulate
try:
    from pydantic.v1 import BaseModel
except ModuleNotFoundError:
    from pydantic import BaseModel
from markdown import Markdown
from markdown.extensions import Extension
from markdown.preprocessors import Preprocessor


class Mdantic(Extension):
    def __init__(self, configs=None):
        if configs is None:
            configs = {}
        self.config = {
            "init_code": ["", "python code to run when initializing"],
            "columns": [
                ["key", "type", "required", "description", "default"],
                "Columns to use in table, comma separated list",
            ],
        }
        for key, value in configs.items():
            self.setConfig(key, value)
        super().__init__()

    def extendMarkdown(self, md: Markdown) -> None:
        md.preprocessors.register(MdanticPreprocessor(md, self.getConfigs()), "mdantic", 100)


Field = namedtuple("Field", "key type required description default")


def analyze(cls_name: str) -> Optional[Dict[str, List[Field]]]:
    paths = cls_name.rsplit(".", 1)
    if len(paths) != 2:
        return None

    module = paths[0]
    attr = paths[1]
    try:
        mod = importlib.import_module(module)
    except ModuleNotFoundError:
        return None
    if not hasattr(mod, attr):
        return None

    cls = getattr(mod, attr)

    if not issubclass(cls, BaseModel):
        return None

    structs = {}
    mk_struct(cls, structs)
    return structs


def get_related_enum(ty: type):
    visited = set()
    result = []

    get_related_enum_helper(ty, visited, result)

    return result


def get_enum_values(e):
    return [x.value for x in list(e)]


def get_related_enum_helper(ty, visited, result):
    visited.add(ty)
    if inspect.isclass(ty) and issubclass(ty, Enum) and ty not in result:
        result.append(ty)

    if hasattr(ty, "__args__"):
        for sub_ty in getattr(ty, "__args__"):
            if sub_ty not in visited:
                get_related_enum_helper(sub_ty, visited, result)


# v1:
def mk_struct(cls: type[BaseModel], structs: Dict[str, List[Field]]) -> None:
    this_struct: List[Field] = []
    structs[cls.__name__] = this_struct
    # v2: for field_name, f in cls.model_fields.items():
    for field_name, f in cls.__fields__.items():
        title = f.field_info.title or field_name
        annotation = str(f.type_)
        description = "" if f.field_info.description is None else f.field_info.description

        if annotation is None:
            return None

        related_enums = get_related_enum(annotation)
        if related_enums:
            for e in related_enums:
                description += f"</br>{e.__name__}: {get_enum_values(e)}"

        default = f.get_default()
        default = None if str(default) == "PydanticUndefined" else str(default)

        if hasattr(annotation, "__origin__"):
            ty = str(annotation)
        elif hasattr(annotation, "__name__"):
            ty = annotation.__name__
        else:
            ty = str(annotation)

        this_struct.append(
            Field(
                title,
                ty,
                # v2: str(f.is_required()),
                str(f.required),
                description,
                default,
            )
        )
        if hasattr(annotation, "__mro__"):
            if BaseModel in annotation.__mro__:
                mk_struct(annotation, structs)


def fmt_tab(structs: Dict[str, List[Field]], columns: List[str]) -> Dict[str, str]:
    tabs = {}
    for cls, struct in structs.items():
        tab = []
        for f in struct:
            tab.append([getattr(f, name) for name in columns])
        tabs[cls] = tabulate.tabulate(tab, headers=columns, tablefmt="github")
    return tabs


class MdanticPreprocessor(Preprocessor):
    """
    This provides an "include" function for Markdown, similar to that found in
    LaTeX (also the C pre-processor and Fortran). The syntax is {!filename!},
    which will be replaced by the contents of filename. Any such statements in
    filename will also be replaced. This replacement is done prior to any other
    Markdown processing. All file-names are evaluated relative to the location
    from which Markdown is being called.
    """

    def __init__(self, md: Markdown, config):
        super(MdanticPreprocessor, self).__init__(md)
        self.init_code = config["init_code"]
        if self.init_code:
            exec(self.init_code)
        self.columns = config["columns"]

    def run(self, lines: List[str]):
        for i, l in enumerate(lines):
            g = re.match(r"^\$pydantic: (.*)$", l)
            if g:
                cls_name = g.group(1)
                structs = analyze(cls_name)
                if structs is None:
                    print(f"warning: mdantic pattern detected but failed to process or import: {cls_name}")
                    continue
                tabs = fmt_tab(structs, self.columns)
                table_str = ""
                for cls, tab in tabs.items():
                    table_str += "\n" + f"**{cls}**" + "\n\n" + str(tab) + "\n"
                lines = lines[:i] + [table_str] + lines[i + 1 :]

        return lines


def makeExtension(*_, **kwargs):
    return Mdantic(kwargs)
