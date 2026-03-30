"""
Some problems/missing features exist in the loading and saving of Functions.externals.MaCSBio_cobrapy models
therefore the necessary functions have been taken from Functions.externals.MaCSBio_cobrapy version 0.29.0
and added to the code.

Author: Jelle Bonthuis, 2025-01-22
Changes made at date:
    - Wrap load_mat method in load_matlab_model to provide useful error message when
        (most likely) non-alpha numeric characters are present in the model
    - Load_matlab_model now loads rxnECNumber as an optional reaction attribute,
        specifically from_mat_struct method is edited to include this attribute.
            EDIT: It appears this is already included in the original code, however
            the m.dtype.names contains eccodes and not rxnECNumber, so this is changed
            to also recognize that.
    - Save_json_model method now deals properly with Na and NaN values in the model by
        setting them to None. However this will provide warnings when doing so,
        indicating that the user should check the model for these values.
    - Save_json_model is now also saving EC_numbers,
        specifically its model_to_dict method.
    - Changed load_mat function to raise errors instead of ignoring them
        when it is impossible to load a specific field.
    - Fixed issue with loading subSystems of reactions from .mat files as this was
        inside a try except block that would ignore the error. The fix included
        not checking for truth value of an array.
    - Added geneShortNames field check, and removed try/except for geneNames
        TODO other saving formats should do so too
"""

# imports
from typing import Union, Optional, IO
from pathlib import Path
import json
import numpy as np
from typing import Dict, Any
from Functions.externals.MaCSBio_cobrapy import Model, Metabolite, Reaction, Gene
from Functions.externals.MaCSBio_cobrapy.io.mat import (
    _cell_to_str_list,
    _cell_to_float_list,
    _get_id_compartment,
    DICT_MET,
    DICT_MET_NOTES,
    DICT_GENE,
    DICT_REACTION,
    DICT_REACTION_NOTES,
    mat_parse_annotations,
    mat_parse_notes,
    set_objective,
    Group,
)
from Functions.externals.MaCSBio_cobrapy.io.dict import (
    _metabolite_to_dict,
    _reaction_to_dict,
    _gene_to_dict,
    _metabolite_from_dict,
    gene_from_dict,
    _OPTIONAL_MODEL_ATTRIBUTES,
    _ORDERED_OPTIONAL_MODEL_KEYS,
    _update_optional,
    _fix_type,
    _reaction_from_dict,
    _REQUIRED_REACTION_ATTRIBUTES,
    attrgetter,
    itemgetter,
)
from Functions.externals.MaCSBio_cobrapy.io.json import JSON_SPEC
from collections import OrderedDict


from logging import getLogger

try:
    from scipy.io.matlab import loadmat
    from scipy.sparse import csc_matrix
    scipy_io_matlab_loadmat = "placeholder_string"
    scipy_sparse_csc_matrix = "placeholder_string"
except ImportError:
    scipy_io_matlab_loadmat = None
    scipy_sparse_csc_matrix = None
logger = getLogger(__name__)

#### Jelle Bonthuis 2025-01-22: Added rxnECNumber to the list of valid fields
_ORDERED_OPTIONAL_REACTION_KEYS = [
    "objective_coefficient",
    "subsystem",
    "notes",
    "annotation",
    # "EC_number",
]
_OPTIONAL_REACTION_ATTRIBUTES = {
    "objective_coefficient": 0,
    "subsystem": "",
    "notes": {},
    "annotation": {},
    # "EC_number": "",
}

_REACTION_FIELDS_TO_CHECK = [
    "id",
    "name",
    "subsystem",
    "gene_reaction_rule",
    "lower_bound",
    "upper_bound",
    "gene_reaction_rule",
]

_META_FIELDS_TO_CHECK = [
    "id",
    "name",
    "formula",
    "charge",
    "compartment",
    "annotation",
    "notes",
    "inchis",
]

_GENE_FIELDS_TO_CHECK = ["id", "name", "annotation", "notes"]


def load_matlab_model(
        infile_path: Union[str, Path, IO],
        variable_name: Optional[str] = None,
        inf: float = np.inf,
) -> Model:
    """Load a Functions.externals.MaCSBio_cobrapy model stored as a .mat file.

    Parameters
    ----------
    infile_path : str or Path or filehandle
        File path or descriptor of the .mat file describing the Functions.externals.MaCSBio_cobrapy model.
    variable_name : str, optional
        The variable name of the model in the .mat file. If None, then the
        first MATLAB variable which looks like a Functions.externals.MaCSBio_cobrapy model will be used
        (default None).
    inf: float, optional
        The value to use for infinite bounds. Some solvers do not handle
        infinite values so for using those, set this to a high numeric value
        (default `numpy.inf`).

    Returns
    -------
    Functions.externals.MaCSBio_cobrapy.Model
        The Functions.externals.MaCSBio_cobrapy model as represented in the .mat file.

    Raises
    ------
    ImportError
        If scipy is not found in the Python environment.
    IOError
        If no Functions.externals.MaCSBio_cobrapy model is found in the .mat file.

    """
    if not scipy_io_matlab_loadmat:
        raise ImportError("load_matlab_model() requires scipy.")
    if not scipy_sparse_csc_matrix:
        raise ImportError("load_matlab_model() requires scipy.")

    ### Jelle Bonthuis 2025-01-22: Try/except added
    try:
        if isinstance(infile_path, str):
            data = loadmat(infile_path)
        elif isinstance(infile_path, Path):
            data = loadmat(infile_path.open("rb"))  # noqa W9018
        else:
            data = loadmat(infile_path)  # noqa W9018
    except TypeError as e:
        if "buffer is too small" in str(e):
            raise RuntimeError(f"Error loading .mat file: {e}"
                              f"This can be due to problems with non-alpha "
                              f"numeric characters (e.g. greek symbols) in the model.")
        else:
            raise IOError(f"Error loading .mat file: {e}")

    possible_names = []
    if variable_name is None:
        # skip meta variables
        meta_vars = {"__globals__", "__header__", "__version__"}
        possible_names = sorted(i for i in data if i not in meta_vars)
        if len(possible_names) == 1:
            variable_name = possible_names[0]
    elif variable_name is not None:
        return from_mat_struct(data[variable_name], model_id=variable_name, inf=inf)

    for possible_name in possible_names:
        # return from_mat_struct(data[possible_name], model_id=possible_name, inf=inf)
        try:
            return from_mat_struct(data[possible_name], model_id=possible_name, inf=inf)
        except ValueError:
            #### added by Jelle Bonthuis 2025-01-22: Added warning for error
            logger.warning(f"Error loading model {possible_name}")
            pass # SHOULD NOT HAVE A RAISE FOR THIS IS ONLY OPTIONAL
            # TODO: use custom Functions.externals.MaCSBio_cobrapy exception to handle exception
    # If code here is executed, then no model was found.
    raise IOError(f"No Functions.externals.MaCSBio_cobrapy model found at {infile_path}.")


def from_mat_struct(
        mat_struct: np.ndarray,
        model_id: Optional[str] = None,
        inf: float = np.inf,
) -> Model:
    """Create a model from the cobratoolbox struct.

    Parameters
    ----------
    mat_struct : numpy.ndarray
        The `numpy.ndarray` that most likely contains the model, being chosen by
        load_matlab_file after loading the matlab structure via scipy_io.loadmat.
    model_id : str, optional
        The ID of the model generated. If None, will try to look for ID in
        model's description. If multiple IDs are found, the first one is
        used. If no IDs are found, will use 'imported_model' (default None).
    inf : float, optional
        The value to use for infinite bounds. Some solvers do not handle
        infinite values so for using those, set this to a high numeric value
        (default `numpy.inf`).

    Returns
    -------
    Functions.externals.MaCSBio_cobrapy.Model
        The model as represented in .mat file.

    """
    m = mat_struct
    #### Jelle Bonthuis 2025-01-22: Added rxnECNumber to the list of valid fields
    ##### as well as eccodes
    if m.dtype.names is None or {"rxns", "mets", "S", "lb", "ub"} > set(m.dtype.names):
        raise ValueError("Invalid MATLAB struct.")

    old_cobratoolbox_fields = [
        "confidenceScores",
        "metCharge",
        "ecNumbers",
        "eccodes",
        "KEGGID",
        "metSmile",
        "metHMDB",
    ]
    new_cobratoolbox_fields = [
        "rxnConfidenceScores",
        "metCharges",
        "rxnECNumber",
        "rxnECNumber",
        "metKEGGID",
        "metSmiles",
        "metHMDBID",
    ]
    for old_field, new_field in zip(old_cobratoolbox_fields, new_cobratoolbox_fields):
        if old_field in m.dtype.names and new_field not in m.dtype.names:
            logger.warning(
                f"This model seems to have {old_field} instead of {new_field} field. "
                f"Will use {old_field} for what {new_field} represents."
            )
            new_names = list(m.dtype.names)
            new_names[new_names.index(old_field)] = new_field
            m.dtype.names = new_names

    model = Model()
    if model_id is not None:
        model.id = model_id
    elif "description" in m.dtype.names:
        description = m["description"][0, 0][0]
        if not isinstance(description, str) and len(description) > 1:
            model.id = description[0]
            logger.warning("Several IDs detected, only using the first.")
        else:
            model.id = description
    else:
        model.id = "imported_model"
    if "modelName" in m.dtype.names and np.size(m["modelName"][0, 0]):
        model.name = m["modelName"][0, 0][0]

    met_ids = _cell_to_str_list(m["mets"][0, 0])
    if {"metComps", "comps", "compNames"}.issubset(m.dtype.names):
        met_comp_index = [x[0] - 1 for x in m["metComps"][0][0]]
        comps = _cell_to_str_list(m["comps"][0, 0])
        comp_names = _cell_to_str_list(m["compNames"][0][0])
        met_comps = [comps[i] for i in met_comp_index]
        met_comp_names = [comp_names[i] for i in met_comp_index]
    else:
        logger.warning(
            f"No defined compartments in model {model.id}. "
            f"Compartments will be deduced heuristically "
            f"using regular expressions."
        )
        met_comps = [_get_id_compartment(x) for x in met_ids]
        met_comp_names = met_comps
        if None in met_comps or "" in met_comps:
            raise ValueError("Some compartments were empty. Check the model!")
        logger.warning(
            f"Using regular expression found the following compartments:"
            f"{', '.join(sorted(set(met_comps)))}"
        )
    if None in met_comps or "" in met_comps:
        raise ValueError("Some compartments were empty. Check the model!")
    model.compartments = dict(zip(met_comps, met_comp_names))
    met_names, met_formulas, met_charges, met_inchis = None, None, None, None
    try:
        met_names = _cell_to_str_list(m["metNames"][0, 0], "")
    except (IndexError, ValueError):
        # TODO: use custom Functions.externals.MaCSBio_cobrapy exception to handle exception
        pass
    try:
        met_formulas = _cell_to_str_list(m["metFormulas"][0, 0])
    except (IndexError, ValueError):
        # TODO: use custom Functions.externals.MaCSBio_cobrapy exception to handle exception
        pass
    try:
        met_charges = _cell_to_float_list(m["metCharges"][0, 0])
    except (IndexError, ValueError):
        # TODO: use custom Functions.externals.MaCSBio_cobrapy exception to handle exception
        pass
    try:
        met_inchis = _cell_to_str_list(m["inchis"][0,0])
    except:
        pass

    new_metabolites = []
    for i in range(len(met_ids)):
        new_metabolite = Metabolite(met_ids[i], compartment=met_comps[i])
        if met_names:
            new_metabolite.name = met_names[i]
        if met_charges:
            new_metabolite.charge = met_charges[i]
        if met_formulas:
            new_metabolite.formula = met_formulas[i]
        if met_inchis:
            new_metabolite.inchi = met_inchis[i]
        new_metabolites.append(new_metabolite)
    mat_parse_annotations(new_metabolites, m, d_replace=DICT_MET)
    mat_parse_notes(new_metabolites, m, d_replace=DICT_MET_NOTES)
    model.add_metabolites(new_metabolites)

    if "genes" in m.dtype.names:
        gene_names = None
        gene_ids = _cell_to_str_list(m["genes"][0, 0])
        ### Jelle Bonthuis 2025-01-22: added geneShortNames field check, and removed try/except
        if "geneNames" in m.dtype.names:
            gene_names = _cell_to_str_list(m["geneNames"][0, 0])
        elif "geneShortNames" in m.dtype.names:
            gene_names = _cell_to_str_list(m["geneShortNames"][0, 0])
        new_genes = [
            Gene(gene_ids[i], name=gene_names[i]) if gene_names else Gene(gene_ids[i])
            for i in range(len(gene_ids))
        ]
        mat_parse_annotations(new_genes, m, d_replace=DICT_GENE)
        for current_gene in new_genes:
            current_gene._model = model
        model.genes += new_genes

    new_reactions = []
    rxn_ids = _cell_to_str_list(m["rxns"][0, 0])
    if "rxnECNumber" in m.dtype.names:
        rxn_EC_number = _cell_to_str_list(m["rxnECNumber"][0, 0], empty_value=None)
    else:
        rxn_EC_number = [None for _ in rxn_ids]
    rxn_lbs = _cell_to_float_list(m["lb"][0, 0], empty_value=None, inf_value=inf)
    rxn_ubs = _cell_to_float_list(m["ub"][0, 0], empty_value=None, inf_value=inf)
    rxn_gene_rules, rxn_names, rxn_subsystems = None, None, None
    try:
        rxn_gene_rules = _cell_to_str_list(m["grRules"][0, 0], "")
    except (IndexError, ValueError) as e:
        ### Jelle Bonthuis 2025-01-22: Added try/except to catch error so that
                # the error will actually show up
        raise e
    try:
        rxn_names = _cell_to_str_list(m["rxnNames"][0, 0], "")
    except (IndexError, ValueError) as e:
        ### Jelle Bonthuis 2025-01-22: Added try/except to catch error so that
        # the error will actually show up
        raise e
    try:
        # RECON3.0 mat has an array within an array for subsystems.
        # If we find a model that has multiple subsytems per reaction, this should be
        # modified
        #### JELLE BONTHUIS 2025-02-11 Added fix as the following try would always fail
                # and no error would be shown, making it so that no subsystems would be loaded

        if np.sctype2char(m["subSystems"][0, 0][0][0]) == "O" and isinstance(
                m["subSystems"][0, 0][0][0][0], np.ndarray
        ):
            ### OLD CODE BEFORE FIX JELE BONTHUIS 2025-02-11
            # rxn_subsystems = [
            #     each_cell[0][0][0][0] if each_cell else ""
            #     for each_cell in m["subSystems"][0, 0]
            # ]

            ### FIX JELE BONTHUIS 2025-02-11
            rxn_subsystems = [
                each_cell[0][0][0][0]
                for each_cell in m["subSystems"][0, 0]
            ]
        # Other matlab files seem normal.
        else:
            rxn_subsystems = _cell_to_str_list(m["subSystems"][0, 0], "")
    except (IndexError, ValueError) as e:
        ### Jelle Bonthuis 2025-01-22: Added try/except to catch error so that
        # the error will actually show up
        raise e
    for i in range(len(rxn_ids)):
        new_reaction = Reaction(
            id=rxn_ids[i],
            lower_bound=rxn_lbs[i],
            upper_bound=rxn_ubs[i],
        )
        if rxn_names:
            new_reaction.name = rxn_names[i]
        if rxn_subsystems:
            new_reaction.subsystem = rxn_subsystems[i]
        if rxn_gene_rules:
            new_reaction.gene_reaction_rule = rxn_gene_rules[i]
        if rxn_EC_number:
            new_reaction.EC_number = rxn_EC_number[i]
        new_reactions.append(new_reaction)
    mat_parse_annotations(new_reactions, m, d_replace=DICT_REACTION)
    # TODO - When cobrapy.notes are revised not be a dictionary (possibly when
    #  annotations are fully SBML compliant, revise this function.
    mat_parse_notes(new_reactions, m, d_replace=DICT_REACTION_NOTES)

    csc = csc_matrix(m["S"][0, 0])
    for i in range(csc.shape[1]):
        stoic_dict = {
            model.metabolites[j]: csc[j, i] for j in csc.getcol(i).nonzero()[0]
        }
        new_reactions[i].add_metabolites(stoic_dict)

    model.add_reactions(new_reactions)

    # Make subsystems into groups
    if rxn_subsystems:
        rxn_group_names = set(rxn_subsystems).difference({None})
        new_groups = []
        for g_name in sorted(rxn_group_names):
            group_members = model.reactions.query(
                lambda x: x.subsystem == g_name  # noqa: B023
            )
            new_group = Group(
                id=g_name, name=g_name, members=group_members, kind="partonomy"
            )
            new_group.annotation["sbo"] = "SBO:0000633"
            new_groups.append(new_group)
        model.add_groups(new_groups)

    if "c" in m.dtype.names:
        c_vec = _cell_to_float_list(m["c"][0, 0])
        coefficients = dict(zip(new_reactions, c_vec))
        if model.solver is not None:
            set_objective(model, coefficients)
    else:
        logger.warning("Objective vector `c` not found.")

    if "osenseStr" in m.dtype.names:
        if isinstance(m["osenseStr"][0, 0][0], np.str_):
            model.objective_direction = str(m["osenseStr"][0, 0][0])
        elif isinstance(m["osenseStr"][0, 0][0], np.ndarray):
            model.objective_direction = str(m["osenseStr"][0, 0][0][0])
    elif "osense" in m.dtype.names:
        osense = float(m["osense"][0, 0][0][0])
        objective_direction_str = "max"
        if osense == 1:
            objective_direction_str = "min"
        model.objective_direction = objective_direction_str
    return model



def _reaction_to_dict(reaction: Reaction) -> OrderedDict:
    """Convert a Functions.externals.MaCSBio_cobrapy Reaction object to a dictionary.

    Parameters
    ----------
    reaction : Functions.externals.MaCSBio_cobrapy.Reaction
        The Functions.externals.MaCSBio_cobrapy.Reaction to convert to dictionary.

    Returns
    -------
    dict
        The converted dictionary object.

    See Also
    --------
    _reaction_from_dict : Convert a dictionary to a Functions.externals.MaCSBio_cobrapy Reaction object.

    """
    new_reaction = OrderedDict()
    for key in _REQUIRED_REACTION_ATTRIBUTES:
        if key != "metabolites":
            if key == "lower_bound" and (
                    np.isnan(reaction.lower_bound) or np.isinf(reaction.lower_bound)
            ):
                new_reaction[key] = str(_fix_type(getattr(reaction, key)))
            elif key == "upper_bound" and (
                    np.isnan(reaction.upper_bound) or np.isinf(reaction.upper_bound)
            ):
                new_reaction[key] = str(_fix_type(getattr(reaction, key)))
            else:
                new_reaction[key] = _fix_type(getattr(reaction, key))
            continue
        mets = OrderedDict()
        for met in sorted(reaction.metabolites, key=attrgetter("id")):
            mets[str(met)] = reaction.metabolites[met]
        new_reaction["metabolites"] = mets
    _update_optional(
        reaction,
        new_reaction,
        _OPTIONAL_REACTION_ATTRIBUTES,
        _ORDERED_OPTIONAL_REACTION_KEYS,
    )
    return new_reaction


def model_to_dict(model: Model, sort: bool = False) -> OrderedDict:
    """Convert a Functions.externals.MaCSBio_cobrapy Model to a dictionary.

    Parameters
    ----------
    model : Functions.externals.MaCSBio_cobrapy.Model
        The model to reformulate as a dict.
    sort : bool, optional
        Whether to sort the metabolites, reactions, and genes or maintain the
        order defined in the model (default False).

    Returns
    -------
    OrderedDict
        A dictionary with keys: 'genes', 'compartments', 'id',
        'metabolites', 'notes' and 'reactions'; where 'metabolites', 'genes'
        and 'metabolites' are in turn lists of dictionaries holding all
        attributes to form the corresponding object.

    See Also
    --------
    model_from_dict : Convert a dictionary to a Functions.externals.MaCSBio_cobrapy Model.

    """
    obj = OrderedDict()
    obj["metabolites"] = list(map(_metabolite_to_dict, model.metabolites))
    obj["reactions"] = list(map(_reaction_to_dict, model.reactions))
    obj["genes"] = list(map(_gene_to_dict, model.genes))
    obj["id"] = model.id
    _update_optional(
        model, obj, _OPTIONAL_MODEL_ATTRIBUTES, _ORDERED_OPTIONAL_MODEL_KEYS
    )
    if sort:
        get_id = itemgetter("id")
        obj["metabolites"].sort(key=get_id)
        obj["reactions"].sort(key=get_id)
        obj["genes"].sort(key=get_id)
    return obj


def model_from_dict(obj: Dict) -> Model:
    """Build a Functions.externals.MaCSBio_cobrapy Model from a dictionary.

    Models stored in JSON are first formulated as a dictionary that can be read
    to a Functions.externals.MaCSBio_cobrapy Model using this function.

    Parameters
    ----------
    obj : dict
        A dictionary with keys: 'genes', 'compartments', 'id',
        'metabolites', 'notes' and 'reactions'; where 'metabolites', 'genes'
        and 'metabolites' are in turn lists of dictionaries holding all
        attributes to form the corresponding object.

    Returns
    -------
    Functions.externals.MaCSBio_cobrapy.Model
        The generated model.

    Raises
    ------
    ValueError
        If `obj` has no 'reactions' attribute.

    See Also
    --------
    model_to_dict : Convert a Functions.externals.MaCSBio_cobrapy Model to a dictionary.

    """
    if "reactions" not in obj:
        raise ValueError("Object has no .reactions attribute. Cannot load.")
    model = Model()
    model.add_metabolites(
        [_metabolite_from_dict(metabolite) for metabolite in obj["metabolites"]]
    )
    model.genes.extend([gene_from_dict(gene) for gene in obj["genes"]])
    model.add_reactions(
        [_reaction_from_dict(reaction, model) for reaction in obj["reactions"]]
    )
    objective_reactions = [
        rxn for rxn in obj["reactions"] if rxn.get("objective_coefficient", 0) != 0
    ]
    coefficients = {
        model.reactions.get_by_id(rxn["id"]): rxn["objective_coefficient"]
        for rxn in objective_reactions
    }
    set_objective(model, coefficients)
    for k, v in obj.items():
        if k in {"id", "name", "notes", "compartments", "annotation"}:
            setattr(model, k, v)
    return model


def save_json_model(
        model: "Model",
        filename: Union[str, Path, IO],
        sort: bool = False,
        pretty: bool = False,
        recursive_call: bool = False,
        **kwargs: Any
) -> None:
    """Write the Functions.externals.MaCSBio_cobrapy model to a file in JSON format.

    Parameters
    ----------
    model : Functions.externals.MaCSBio_cobrapy.Model
        The Functions.externals.MaCSBio_cobrapy model to represent.
    filename : str or file-like
        File path or descriptor that the JSON representation should be
        written to.
    sort : bool, optional
        Whether to sort the metabolites, reactions, and genes or maintain the
        order defined in the model (default False).
    pretty : bool, optional
        Whether to format the JSON more compactly (default) or in a more
        verbose but easier to read fashion. Can be partially overwritten by the
        `**kwargs` (default False).
    **kwargs : Any
        Keyword arguments passed to `json.dump`.

    See Also
    --------
    to_json : Return a string representation.
    json.dump : Base function.

    """
    obj = model_to_dict(model, sort=sort)
    obj["version"] = JSON_SPEC

    if pretty:
        dump_opts = {
            "indent": 4,
            "separators": (",", ": "),
            "sort_keys": True,
            "allow_nan": False,
        }
    else:
        dump_opts = {
            "indent": 0,
            "separators": (",", ":"),
            "sort_keys": False,
            "allow_nan": False,
        }
    dump_opts.update(**kwargs)

    #### Jelle Bonthuis 2025-01-22: Added try/except to catch Na and NaN values
    try:
        if isinstance(filename, (str, Path)):
            with open(filename, "w") as file_handle:
                json.dump(obj, file_handle, **dump_opts)
        else:
            json.dump(obj, filename, **dump_opts)
    except ValueError as e:
        if "Out of range float values are not JSON compliant: nan" in str(e):
            if recursive_call:
                raise RuntimeError("Error saving model to JSON: NaN values in model."
                                   "And already tried to fix them. But still not JSON compliant.")
            logger.warning(
                "Some values in the model are not JSON compliant. "
                "These values will be set to None. "
                "Please check the model for these values."
            )
            for reaction in model.reactions:
                # do this for every field of the reaction
                for field in reaction.__dict__:
                    try:
                        if (
                            field in _REACTION_FIELDS_TO_CHECK
                            and not isinstance(getattr(reaction, field), str)
                            and (getattr(reaction, field) in ["Na", "NaN"]
                            or np.isnan(getattr(reaction, field)))
                        ):
                            setattr(reaction, field, None)
                    except:
                        pass

            for metabolite in model.metabolites:
                for field in metabolite.__dict__:
                    try:
                        if (
                            field in _META_FIELDS_TO_CHECK
                            and not isinstance(getattr(metabolite, field), str)
                            and (getattr(metabolite, field) in ["Na", "NaN"]
                            or np.isnan(getattr(metabolite, field)))
                        ):
                            setattr(metabolite, field, None)
                    except:
                        pass

            for gene in model.genes:
                for field in gene.__dict__:
                    try:
                        if (
                            field in _GENE_FIELDS_TO_CHECK
                            and not isinstance(getattr(gene, field), str)
                            and (getattr(gene, field) in ["Na", "NaN"]
                            or np.isnan(getattr(gene, field)))
                        ):
                            setattr(gene, field, None)
                    except:
                        pass
            save_json_model(model,
                            filename,
                            sort=sort,
                            pretty=pretty,
                            recursive_call=True,
                            **kwargs)
        else:
            raise RuntimeError(f"Error saving model to JSON: {e}")

        if "of type int64 is not JSON serializable" in str(e):
            if recursive_call:
                raise RuntimeError("Error saving model to JSON: int64 values in model."
                                   "And already tried to fix them. But still not JSON compliant.")
            logger.warning(
                "Some values in the model are not JSON (int64).\n"
                "These values will be set to int as json is not numpy compatible.\n"
                "This might take a while"
            )
            for reaction in model.reactions:
                # do this for every field of the reaction
                for field in reaction.__dict__:
                    try:
                        if (
                            field in _REACTION_FIELDS_TO_CHECK
                            and not isinstance(getattr(reaction, field), str)
                            and isinstance(getattr(reaction, field), np.int64)
                        ):
                            setattr(reaction, field, int(getattr(reaction, field)))
                    except:
                        pass
            for metabolite in model.metabolites:
                for field in metabolite.__dict__:
                    try:
                        if (
                            field in _META_FIELDS_TO_CHECK
                            and not isinstance(getattr(metabolite, field), str)
                            and isinstance(getattr(metabolite, field), np.int64)
                        ):
                            setattr(metabolite, field, int(getattr(metabolite, field)))
                    except:
                        pass
            for gene in model.genes:
                for field in gene.__dict__:
                    try:
                        if (
                            field in _GENE_FIELDS_TO_CHECK
                            and not isinstance(getattr(gene, field), str)
                            and isinstance(getattr(gene, field), np.int64)
                        ):
                            setattr(gene, field, int(getattr(gene, field)))
                    except:
                        pass
        else:
            raise RuntimeError("Error saving model to JSON: int64 values in model."
                               "These values are not JSON serializable. "
                               "Please check the model for these values.")

