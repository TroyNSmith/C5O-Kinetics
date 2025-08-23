import autode as ade

ade.Config.n_cores = 16

print("n_cores:", ade.Config.n_cores)
print("get_lmethod:", ade.methods.get_lmethod())
print("get_hmethod:", ade.methods.get_hmethod())

rxn = ade.Reaction("C.[OH]>>[CH3].O", name="test")

rxn.calculate_reaction_profile()
