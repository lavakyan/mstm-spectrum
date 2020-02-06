from mstm_studio.alloy_AuAg import AlloyAuAg

au1ag2 = AlloyAuAg(x_Au=1./3)
fig, axs = au1ag2.plot()
fig.savefig('mat_au1ag2.png', bbox_inches='tight')
