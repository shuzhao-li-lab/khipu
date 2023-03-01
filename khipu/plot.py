try:
    import matplotlib.pyplot as plt
except:
    print("[khipu.utils] cannot import matplotlib, skipping.")


def plot_khipugram(df, savepdf='', relabel={}):
    '''Plot the khipu grid as diagram.
    df = KP.get_khipu_intensities()
    The trunk labels are updated here but not in KP instances, because
    we use the search adduct tables to index m/z distances to M+H+,
    and it's easier to keep the names consistent with the input tables.  
    '''
    _M, _N = df.shape
    zdata = []
    for ii in range(_M):
        for jj in range(_N):
            zdata.append((jj, ii, df.iloc[ii, jj]))

    X = [d[0] for d in zdata]
    Y = [d[1] for d in zdata]
    S = [(np.log10(d[2]+1))**2 for d in zdata]

    trunk_labels = list(df.columns)
    _base = trunk_labels[0]
    trunk_labels = [_base] + [_base+","+x for x in trunk_labels[1:]]
    if relabel:
        trunk_labels = [relabel.get(x, x) for x in trunk_labels]

    fig, ax = plt.subplots()
    for jj in range(_N):
        ax.text(jj, -1, trunk_labels[jj], rotation=60)
        ax.plot([jj]*_M, range(_M), marker='o', linestyle='--', markersize=0.1)

    ax.plot([-1, _N+1], [0,0], linestyle='-', linewidth=2, color='k', alpha=0.3)
    ax.scatter(X, Y, c='red', s=S, alpha=0.8)
    for ii in range(_M):
        ax.text(_N+1.6, ii, df.index[ii])

    ax.margins(0.2)
    ax.set_axis_off()
    ax.invert_yaxis()
    fig.tight_layout()
    if savepdf:
        plt.savefig(savepdf)
    else:
        plt.show()
