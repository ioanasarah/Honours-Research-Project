import cobra 
from cobra import Model, Reaction, Metabolite
from cobra.io import (load_matlab_model)


matlab_file = r"C:\Ioana\_uni\honours\slp\algorithm\Main_files\models\updated_paracetamol_specific\model_Human_GEM_with_paracetamol_reactions.mat"
model = load_matlab_model(matlab_file)
print('loaded model')

def adding_reactions():

    # create a new reaction
    reaction = cobra.Reaction('MAR30000')
    reaction.name = 'Detoxification of NAPQI using GSH'
    # reaction.subsystem = ''
    reaction.lower_bound = 0.0
    reaction.upper_bound = 1000.0
    reaction.gene_reaction_rule = '(ENSG00000134184 or ENSG00000084207 or ENSG00000277656)'
    print('created reaction')

    reaction_cys = cobra.Reaction('MAR30001')
    reaction_cys.name = 'APAP-GSH to APAP-Cys'
    reaction_cys.lower_bound = 0.0
    reaction_cys.upper_bound = 1000.0
    reaction_cys.gene_reaction_rule = '(ENSG00000100031 or ENSG00000197635)'
    print('created reaction_cys')

    reaction_mer = cobra.Reaction('MAR30002')
    reaction_mer.name = 'APAP-GSH to APAP-Mer'
    reaction_mer.lower_bound = 0.0
    reaction_mer.upper_bound = 1000.0
    reaction_mer.gene_reaction_rule = '(ENSG00000171428 or ENSG00000156006)'
    print('created reaction_mer')

    reaction_cys_to_mer = cobra.Reaction('MAR30003')
    reaction_cys_to_mer.name = 'APAP-Cys to APAP-Mer'
    reaction_cys_to_mer.lower_bound = 0.0
    reaction_cys_to_mer.upper_bound = 1000.0
    reaction_cys_to_mer.gene_reaction_rule = '(ENSG00000171428 or ENSG00000156006)'
    print('created reaction_cys_to_mer')


    # creating APAP-GSH metabolite
    APAPGSH = Metabolite(
    'MAM10000c',
    formula='C18H24N4O8S',
    name='Acetaminophen-glutathione',
    compartment='c')
    print('created APAPGSH metabolite')

    # creating mercapturic acid metabolite
    mercapturate = Metabolite(
    'MAM10001c',
    formula='C5H9NO3S',
    name='Mercapturic acid',
    compartment='c')
    print('created mercapturate metabolite')


    # add metabolites
    reaction.add_metabolites({
        model.metabolites.get_by_id("MAM03779r"): -1.0,  # N-acetyl-p-benzoquinone imine
        model.metabolites.get_by_id('MAM02026r'): -1.0,    # Glutathione
        APAPGSH: 1.0,   
    })

    reaction_cys.add_metabolites({
        APAPGSH: -1.0,  # Acetaminophen-glutathione
        model.metabolites.get_by_id('MAM01628c'): -1.0,    # Cysteine
        # model.metabolites.get_by_id('MAM02553c'): -2.0,    # NADPH
        model.metabolites.get_by_id('MAM03526c'): 1.0,    # APAP-cysteine
        model.metabolites.get_by_id('MAM02026c'): 1.0    # GSH
        # model.metabolites.get_by_id('MAM02554c'): 2.0,    # NADP+
        # model.metabolites.get_by_id('MAM02040i'): 1.0,    # h2o
    })

    reaction_mer.add_metabolites({
        APAPGSH: -1.0,    # APAP-GSH  
        mercapturate: -1.0,    # meracpturate
        model.metabolites.get_by_id('MAM03759c'): 1.0,    # APAP-mer
        model.metabolites.get_by_id('MAM02026c'): 1.0,    # GSH
    })
    reaction_cys_to_mer.add_metabolites({
        model.metabolites.get_by_id('MAM03526c'): -1.0,    # APAP-Cys
        mercapturate: -1.0,    # meracpturate
        model.metabolites.get_by_id('MAM03759c'): 1.0,    # APAP-mer
        model.metabolites.get_by_id('MAM01628c'): 1.0,    # cysteine

    })
    return reaction, reaction_cys, reaction_mer, reaction_cys_to_mer

def adding_gene_rules():
    reaction_napqi_orig = model.reactions.get_by_id('MAR13001')
    reaction_napqi_orig.gene_reaction_rule = '(ENSG00000140505)' # CYP2A1

def adding_correct_napqi_formation():
    reaction_napqi  = cobra.Reaction('MAR30004')
    reaction_napqi.name = 'Formation of NAPQI directly from cytosolic Acetaminophen'
    reaction_napqi.lower_bound = 0.0
    reaction_napqi.upper_bound = 1000.0
   

    reaction_napqi.add_metabolites({
        model.metabolites.get_by_id('MAM03402c'): -1.0,    # APAP[c]
        model.metabolites.get_by_id('MAM02555c'): -1.0,    # NADPh[c]
        model.metabolites.get_by_id('MAM02630c'): -1.0,    # O2[c]
        model.metabolites.get_by_id('MAM02039c'): -1.0,    # H+[c]
        model.metabolites.get_by_id('MAM03779r'): 1.0,  # N-acetyl-p-benzoquinone imine [r]
        model.metabolites.get_by_id('MAM02554r'): 1.0,    # NADP+[r]
        model.metabolites.get_by_id('MAM02040r'): 2.0,    # H2O[r]
    })
    reaction_napqi.gene_reaction_rule = '(ENSG00000130649 or ENSG00000160868)' # CYP2E1(major) or CYP3A4(minor)
    print('created reaction_napqi')
    print(reaction_napqi)
    return reaction_napqi

def updating_glucuronidation():
    reaction_glucurinidation  = cobra.Reaction('MAR30005')
    reaction_glucurinidation.name = 'Glucoronidation of Acetaminophen (using memberane proteins)'
    reaction_glucurinidation.lower_bound = 0.0
    reaction_glucurinidation.upper_bound = 1000.0
   

    reaction_glucurinidation.add_metabolites({
        model.metabolites.get_by_id('MAM03402c'): -1.0,    # APAP[c]
        model.metabolites.get_by_id('MAM03109c'): -1.0,    # UDP-glucuronate [c]
        model.metabolites.get_by_id('MAM03106r'): 1.0,     #UDP [r]
        model.metabolites.get_by_id('MAM02039r'): 1.0,    # H+[r]
        model.metabolites.get_by_id('MAM03403r'): 1.0,    #Acetaminophen glucuronide
    })
    reaction_glucurinidation.gene_reaction_rule = '(ENSG00000241635 or ENSG00000167165 or ENSG00000241119 or ENSG00000196620)' # UGT1A1 or UGT1A6 or UGT1A9 or UGT2B15 
    print('created reaction_glucurinidation')
    print(reaction_glucurinidation)
    return reaction_glucurinidation


if __name__ == "__main__":
    # reaction, reaction_cys, reaction_mer, reaction_cys_to_mer = adding_reactions()
    # model.add_reactions([reaction, reaction_cys, reaction_mer, reaction_cys_to_mer])
    # reaction_napqi = adding_correct_napqi_formation()
    # model.add_reactions([reaction_napqi])
    # model.reactions.get_by_id('MAR13001').remove_from_model() 
    model.reactions.get_by_id('MAR30005').remove_from_model()
    print('removed reactions')
    reaction_glucurinidation = updating_glucuronidation()
    model.add_reactions([reaction_glucurinidation])
    print('added updated glucuronidation reaction')
    # model.reactions.get_by_id('MAR12462').remove_from_model()
    # print('added reactions to model')
    # print(reaction)
    # print(reaction_cys)
    # print(reaction_mer)
    # print(reaction_cys_to_mer)
    # save the model

    cobra.io.save_matlab_model(model,r"C:\Ioana\_uni\honours\slp\algorithm\Main_files\models\updated_paracetamol_specific\model_Human_GEM_with_paracetamol_reactions.mat")
    print('saved model')