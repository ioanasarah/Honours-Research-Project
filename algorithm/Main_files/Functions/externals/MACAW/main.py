# main.py
'''
TAKEN FROM https://github.com/Devlin-Moyer/macaw version 1.0.0
(commit 7aa9c01e071b3ffb4ab5a8b8cb80db7a70af05f2) as this appears not to be pip-installable
so instead I (Jelle Bonthuis) added it manually to the externals folder on 2025-27-02
After a PR was made on 2025-03-04, a fix was added to the utils.py which has been
manually added here (commit 667347a0c4454c088240399079d824317dec2c3a)

Used in and created for https://www.biorxiv.org/content/10.1101/2024.06.24.600481v1.full
Moyer, Devlin C., et al. "Semi-Automatic Detection of Errors in Genome-Scale Metabolic Models." bioRxiv (2024): 2024-06.

Perform a series of tests on a given Cobrapy Model object to identify reactions
that may have some errors in their definition/construction (e.g. duplicates of
other reactions, incorrect reversibility, part of dead-ends) and suggest ways to
fix at least some of these (potential) errors
'''

from Functions.externals.MACAW.structural import dead_end_test, duplicate_test, diphosphate_test
from Functions.externals.MACAW.flux_based import loop_test, dilution_test
from Functions.externals.MACAW.utils import simplify_test_results, add_reaction_equations
from Functions.externals.MaCSBio_cobrapy.util import Zero
import os
import json
import networkx as nx
import pandas as pd

def run_all_tests(
    model,
    # only required for the redox portion of duplicate_test; will just skip that
    # bit without raising an exception if not provided
    redox_pairs = list(), proton_ids = list(),
    # only required for the reversible diphosphate test; will just skip without
    # raising an exception if not provided
    diphosphate_met_ids = list(), phosphate_met_ids = list(),
    # optional parameters for dilution_test and/or loop_test
    media_mets = None, timeout = 1800, max_attempts = 3, zero_thresh = 10**-8,
    corr_thresh = 0.9,
    # arguments affecting how the reaction equation column in the table of test
    # results are created: with metabolite IDs (default) or names (use_names =
    # True), and whether or not to add suffixes to those IDs or names indicating
    # the compartment the metabolite is in (default is to not)
    use_names = False, add_suffixes = False,
    # only the dilution and loop tests can be run in parallel, but they can
    # really benefit from multiple threads
    threads = 1,
    # control whether or not the tests print progress updates (0 to silence)
    verbose = 1,
    save_location: str|None = None
):
    # identify all dead-end reactions and reactions that are nominally
    # reversible but only actually capable of sustaining flux in a single
    # direction
    (dead_end_results, dead_end_edges) = dead_end_test(
        model, use_names, add_suffixes, verbose
    )
    # print row 6 of dead_end_results
    print(dead_end_results.iloc[6])

    if save_location is not None:
        dead_end_results.to_csv(os.path.join(save_location, "dead_end_results.csv"))
        with open(os.path.join(save_location, "dead_end_edges.json"), 'w') as f:
            json.dump(dead_end_edges, f)

    # # IMPORTANT REMOVE
    # # TODO REMOVE
    # # maximize flux through MAR03905, MAR03907, MAR05000, MAR09099, MAR08756, MAR09242
    met = model.metabolites.get_by_id("MAM02403x")
    print(met.reactions)
    # ids = ["3007 4281 8977 "]
    ids = ["MAR03007", "MAR04281"]
    rxns = [model.reactions.get_by_id(id) for id in ids]
    # max flux through 3007
    model.objective = Zero
    rxns[0].objective_coefficient = 1
    model.objective.direction = "max"
    model.optimize()



    print("placeholder")
    # reactions_to_maximize = ["MAR03907"] #"MAR05000", "MAR09099", "MAR08756", "MAR09242","MAR03905", "MAR04948"]
    # reversibility_of_reactions = {reaction_id: model.reactions.get_by_id(reaction_id).reversibility for reaction_id in reactions_to_maximize}
    # # set objective to 0
    # model.objective = Zero
    # # set objective to maximize flux through each reaction
    # for reaction_id in reactions_to_maximize:
    #     model.reactions.get_by_id(reaction_id).objective_coefficient = 1
    # model.objective.direction = "max"
    # model.optimize()
    #
    # reactions_to_check = ["MAR03905", "MAR03907", "MAR05000", "MAR09099", "MAR08756", "MAR09242","MAR03905", "MAR04948"]
    # # get fluxes for each reaction
    # any_non_zero_fluxes = {reaction_id: model.reactions.get_by_id(reaction_id.id).flux for reaction_id in model.reactions if model.reactions.get_by_id(reaction_id.id).flux != 0}
    # fluxes = {reaction_id: model.reactions.get_by_id(reaction_id).flux for reaction_id in reactions_to_check}
    # names = {reaction_id: model.reactions.get_by_id(reaction_id).name for reaction_id in reactions_to_check}
    # reaction_formulae = {reaction_id: model.reactions.get_by_id(reaction_id).build_reaction_string(use_metabolite_names = True) for reaction_id in reactions_to_check}
    # metabolites_in_reaction = {reaction_id: model.reactions.get_by_id(reaction_id).metabolites for reaction_id in reactions_to_check}
    # metabolite_objects = {reaction_id: {metabolite_id: model.metabolites.get_by_id(metabolite_id.id) for metabolite_id in model.reactions.get_by_id(reaction_id).metabolites} for reaction_id in reactions_to_check}
    #
    # # print reaction names with flux
    # for reaction_id in reactions_to_check:
    #     print("\n")
    #     print(f"{reaction_id}, Flux: {fluxes[reaction_id]} Formula: {reaction_formulae[reaction_id]}")
    #     print(f"Metabolites: {metabolites_in_reaction[reaction_id]}")
    # # TODO END OF REMOVE



    # identify sets of reactions that are potentially duplicates of each other
    (duplicates, dupe_edges) = duplicate_test(
        model, redox_pairs, proton_ids, use_names, add_suffixes, verbose
    )
    if save_location is not None:
        duplicates.to_csv(os.path.join(save_location, "duplicates.csv"))
        with open(os.path.join(save_location, "dupe_edges.json"), 'w') as f:
            json.dump(dupe_edges, f)
    # identify all reversible reactions that involve diphosphate that aren't
    # transporting it between compartments cuz they should probably be
    # irreversible (only if lists of IDs of metabolites representing diphosphate
    # and phosphate were provided)
    diphosphates = diphosphate_test(
        model, diphosphate_met_ids, phosphate_met_ids, use_names, add_suffixes,
        verbose
    )
    if save_location is not None:
        diphosphates.to_csv(os.path.join(save_location, "diphosphates.csv"))

    # identify reactions that are capable of sustaining non-zero fluxes when all
    # exchange reactions are blocked
    (loops, loop_edges) = loop_test(
        model, zero_thresh, corr_thresh, use_names, add_suffixes, threads,
        verbose
    )
    if save_location is not None:
        loops.to_csv(os.path.join(save_location, "loops.csv"))
        with open(os.path.join(save_location, "loop_edges.json"), 'w') as f:
            json.dump(loop_edges, f)
    # identify reactions that become incapable of sustaining non-zero fluxes
    # when dilution constraints are added to the model
    (dilution_results, dilution_edges) = dilution_test(
        model, dead_end_results, media_mets, zero_thresh = zero_thresh,
        timeout = timeout, max_attempts = max_attempts, use_names = use_names,
        add_suffixes = add_suffixes, verbose = verbose, threads = threads,
        save_location = save_location
    )
    # merge the dataframes containing the results of all tests (dead-end and
    # dilution test results were already merged into dilution_test_results)
    all_test_results = duplicates.merge(diphosphates).merge(loops)
    all_test_results = all_test_results.merge(dilution_results)
    # if verbose isn't 0, print the number of reactions flagged by at least one
    # test
    simple_results = simplify_test_results(all_test_results)
    flagged_rxns = simple_results[
        simple_results.loc[
            :, simple_results.columns.str.contains('test')
        ].apply(lambda col: col != 'ok', axis = 1).any(axis = 1)
    ]['reaction_id']
    if verbose > 0:
        msg = f'{len(flagged_rxns)} of the {len(model.reactions)} reactions in '
        msg += 'the given GSMM were flagged by at least one of the tests.'
        print(msg)
    # combine the edge lists and add edges between reactions flagged by
    # different tests. also add a column to the test results indicating which
    # "pathway" (connected component of the network) each reaction winds up in
    (all_test_results, edge_list) = form_pathways(
        all_test_results, model, dilution_edges + dead_end_edges,
        dupe_edges + loop_edges, use_names, add_suffixes
    )
    return((all_test_results, edge_list))

def form_pathways(
    test_results, model, met_rxn_edges, rxn_rxn_edges,
    use_names = False, add_suffixes = False
):
    '''
    Given:
        - a set of tuples of metabolite and reaction IDs generated by the
          dead-end and/or dilution tests (or an empty list)
        - a set of tuples of pairs of reaction IDs generated by the duplicate
          test (or an empty list)
        - the results of one or more other tests
        - the Cobrapy Model object those tests were run on
    Return a Pandas Dataframe with a "source" and "target" column representing
    an edge list of a network connecting reactions flagged by the various tests.
    Note that many reactions flagged by one or more tests may not wind up
    connected to any other reactions.
    '''
    # start by getting all metabolite and reaction IDs that appear in the given
    # metabolite-reaction edge list
    mets = {e[0] for e in met_rxn_edges}
    rxns = {e[1] for e in met_rxn_edges}
    # then add new edges to the metabolite-reaction edge list between all
    # reactions that were flagged by any test and metabolites they involve that
    # are already present in that edge list
    simple_results = simplify_test_results(test_results)
    flagged_reactions = simple_results[simple_results.loc[
        :, simple_results.columns.str.contains('test')
    ].apply(lambda col: col != 'ok', axis = 1).any(axis = 1)]['reaction_id']
    rxns.update(flagged_reactions.to_list())
    # also add edges to any reactions in the reaction-reaction edge lists
    rxns.update({e[0] for e in rxn_rxn_edges})
    rxns.update({e[1] for e in rxn_rxn_edges})
    # coerce to set so we can avoid adding duplicate edges easily
    met_rxn_edges = set(met_rxn_edges)
    for rxn_id in rxns:
        for met in model.reactions.get_by_id(rxn_id).metabolites:
            if (met.id in mets):
                met_rxn_edges.add((met.id, rxn_id))
    # replace all the metabolite-reaction edges with reaction-reaction edges
    # between all reactions that share at least one metabolite (only considering
    # the metabolites present in the metabolite-reaction edge list to avoid
    # connecting all reactions that involve water, protons, ATP, etc.)
    bip_graph = nx.from_edgelist(met_rxn_edges)
    rxn_graph = nx.algorithms.bipartite.projected_graph(
        bip_graph, [n for n in bip_graph if n in rxns]
    )
    edge_df = nx.to_pandas_edgelist(rxn_graph)
    # now we can concatenate this edge list with the given reaction-reaction
    # edge list
    edge_df = pd.concat([
        edge_df, nx.to_pandas_edgelist(nx.from_edgelist(rxn_rxn_edges))
    ]).drop_duplicates()
    # enumerate the connected components of this graph and add a column to the
    # given test results with these numbers so that reactions can be grouped
    # by the "pathway" that they're in
    pathway_dict = {r.id : 0 for r in model.reactions}
    i = 0
    for pathway in nx.connected_components(rxn_graph):
        i += 1
        for reaction_id in pathway:
            pathway_dict[reaction_id] = i
    test_results['pathway'] = test_results['reaction_id'].map(pathway_dict)
    return((test_results, edge_df))
