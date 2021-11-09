def plot_mut_spectrum(freqs, seq_amino, title, ref = "TEM", data = "freq"): 
    #to plot MSA mutational spectrum
    #-------USEFUL
    
    # fix the sides
    pos1 = 1
    pos2 = len(seq_amino)
        
    # length of the sequence to be displayed
    l = pos2 - pos1 + 1
    
    # sequence to display
    seq_amino = seq_amino[(pos1 -1): (pos2)]
    
    # sequence with None instead of gaps
    seq_amino_none_gap = [None if val == 20 else val for val in seq_amino]
    
    # add None in the positions where the wt has a gap
    freqs_nan = np.copy(freqs)
    for i, val in enumerate(seq_amino_none_gap):
        if val == None:
            freqs_nan[:, i] = [None for i in range(20)]

    #-------COLORS
    
    # define colors
    if data == "freq":
        norm1 = matplotlib.colors.LogNorm(vmin = 10^(-3), vmax = 1)
        labb = "Frequency" 
        ccc = ["white",  "#6a72da", "navy",  "limegreen", "limegreen"]
        cvals = [0, 0.16, 0.3, 0.5, 1] 
    elif data == "en":
        ccc = ["black", "green", "limegreen", "white","white", "#6a72da", "darkblue"]
        norm1 = norm1 = matplotlib.colors.Normalize(-6, 10)
        labb = "ΔE"
        cvals = map(norm1, [-6, -5, -1,  1.5, 4.5 , 7, 10])     
    elif data == "fit":
        ccc = ["darkblue", "#2932a6", "white", "limegreen",  "green"]
        norm1 = matplotlib.colors.LogNorm( 0.01, 2.5 )
        labb = "Fitness"
        cvals = map(norm1, [0.01, 0.1, 1.6, 2, 2.5]) 
    
    tuples =list(zip(cvals, ccc))
    cmap1 = matplotlib.colors.LinearSegmentedColormap.from_list("", tuples, N = 50)
        
    # useful stuff
    bar_dict = dict([ ("shrink", 1.2), ("orientation", "vertical")])

    #-------MAINFIG
    
    # start figure
    fig1 = plt.figure(figsize=(25,3), dpi = 300)
    dat =  freqs_nan[:, range(pos1-1, pos2)]
    ax = sns.heatmap(dat, cmap = cmap1,  norm = norm1,
        square = True, center = 0.9,  linewidths=0.01, linecolor = "lightgrey", 
        cbar= True, cbar_kws = bar_dict)
    # add title
    plt.title(title, fontsize=20, y=1.08)
    
    # add labels to axis
    plt.ylabel("Amino acid",rotation=90,fontsize=20)
    plt.xlabel("Protein sequence position",rotation=0,fontsize=20) 

    # make frame visible
    for key in ["top", "bottom", "left", "right"]:
        ax.spines[key].set_visible(True)
        
    #-------COLORBAR
    
    # colorbal labels
    cbar = ax.collections[0].colorbar
    if data != "en":
        cbar.set_ticks([0, 10**(-3), 10**(-2), 10**(-1),0.5, 1, 2])
        cbar.set_ticklabels([0, "10⁻³", "10⁻²", "10⁻¹", 0.5, 1, 2])      
    else:
        cbar.set_ticks(range(-10, 10 + 1, 2))
        cbar.set_ticklabels(range(-10, 10 + 1, 2)) 
    cbar.set_label(labb,fontsize = 18)
    cbar.ax.tick_params(labelsize=18)
    ax.set_facecolor("grey")
    cbar.outline.set_visible(True)
    cbar.outline.set_edgecolor("black")
    cbar.outline.set_linewidth(1)
    axbar = cbar.ax
    axbar.set_aspect(13)
    
    #-------WT
    
    # add red scatter for wt
    amino_pos = [None if amino == None else amino + 0.5 for amino in seq_amino_none_gap]
    ok_pos = [None if val == None else i + 0.5 for i, val in enumerate(seq_amino_none_gap)]
    p = plt.scatter(ok_pos, amino_pos, color = "red", s = 5, alpha = 0.8, 
        linewidth = 0.0, edgecolors = "red")
    
    # add wt name  
    pos_leg = 1.35
    
    if ref == "AAC6":
        pos_leg = 1.25
    lgn = p.axes.legend([f"{ref}"] ,
        bbox_to_anchor = (1,pos_leg), 
    frameon = True, framealpha = 0, edgecolor = "black", fontsize = 15)
    lgn.legendHandles[0]._sizes = [130]
  
        
     #-------XLABELS
    
    # xlabels and xticks
    xxx_pos = [i for i in np.arange( 0.5, l - 0.5) ] 
    labs = [i for i in range(pos1, pos2+1)]
    ax.set_xticks(xxx_pos)
    ax.set_xticklabels(labels = labs, rotation=0, horizontalalignment="right", 
        color="k",fontsize= 12, ha = "center")
    
    # set xtiks every 10
    temp = ax.xaxis.get_major_ticks()
    for i, label in enumerate(temp): 
        if (i + 10)%10 != 0: 
            label.set_visible(False)
            
            
    #-------YLABELS
    
    # set ylabels to nothing
    ax.set_yticks([] )

    # add ylabel
    pos_amino_letts = 208
    amino_dists = 1.8
    amino_start = -8
    if ref == "AAC6":
        amino_dists = 1.2
    if ref == "AAC6":
        pos_amino_letts = 121
    if ref == "AAC6":
        amino_start = -2
    yyy = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L",  "M",  "N", "P",  "Q",  "R",
    "S",  "T", "V",  "W",  "Y"]
    for i in range(20):
        pp = amino_start + amino_dists*i + 1.5
        ax.text(pos_amino_letts, pp, yyy[i], fontsize = 12)
        
    # draw one black bar
    ax2 = plt.axes([0.745, 0.75, 0.018, 0.2])
    x_values = [ 0, 1]
    y_values = [0, 3]
    ax2.plot(x_values, y_values, color = "black")
    ax2.axis("off")
    
    # draw the other black bar
    ax3 = plt.axes([0.745, 0.05, 0.018, 0.2])
    x_values = [ 0, 1]
    y_values = [0, -3]
    ax3.plot(x_values, y_values, color = "black")
    ax3.axis("off")
    return fig1
