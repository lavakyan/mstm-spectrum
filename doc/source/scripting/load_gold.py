from mstm_studio.mstm_spectrum import Material

gold = Material('etaGold.txt')
fig, axs = gold.plot()
fig.savefig('loaded_gold.png', bbox_inches='tight')
