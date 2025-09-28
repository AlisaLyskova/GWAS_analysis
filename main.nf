params.vcf = 'test_data.vcf'
params.phenotype = 'test_data.tsv'
params.output_dir = 'results'

process GWAS_analysis {
    input:
    file(vcf)
    file(phenotype)
    path output_dir

    output:
    path "${output_dir}/*"

    script:
    """
    #!/usr/bin/env python3
    import sys
    import pandas as pd
    import numpy as np
    import cyvcf2
    import allel
    import os
    import matplotlib.pyplot as plt
    import seaborn as sns

    from sklearn.decomposition import PCA
    from scipy.spatial.distance import pdist, squareform
    from scipy.cluster import hierarchy
    from scipy import stats
    from natsort import natsorted
    from matplotlib.colors import LinearSegmentedColormap

    # Functions
    def hamming_distance(u, v):
        #excluding unknown genotypes
        mask = (u != 3) & (v != 3)
        if np.sum(mask) == 0:
            return 1.0
        return np.sum(u[mask] != v[mask]) / np.sum(mask)

    def get_xticks_and_labels(chrom_pos, chrom_sizes):
        ticks = []
        labels = []
        for chrom in chrom_pos:
            start = chrom_pos[chrom]
            end = start + chrom_sizes[chrom]
            ticks.append((start + end) / 2)
            labels.append(chrom)
        return ticks, labels


    # INPUT files
    vcf_file = '${vcf}'
    pheno_file = '${phenotype}'
    output_dir = '${output_dir}'

    colors = ['#FF6700', '#36413E', '#47A025', '#C0C0C0', '#6B0504']
    accent_colors = ['#FF6700', '#6B0504', '#47A025', '#36413E']
    base_color = "#C0C0C0"
    cmap = LinearSegmentedColormap.from_list('custom_gradient', colors)


    pheno = pd.read_csv(pheno_file, sep="\t")
    pheno = pheno.rename(columns={'sample': 'samples'})
    pheno["phenotype"] = pheno["phenotype"].astype("str")

    # Reading vcf
    vcf = cyvcf2.VCF(vcf_file)
    samples = vcf.samples
    genotypes_samples = []
    for record in vcf:
        genotypes_array = list(record.gt_types)
        genotypes_samples.append(genotypes_array)
    df_genotype_by_sample = pd.DataFrame(genotypes_samples, columns=samples)


    # PCA
    pca = PCA(n_components=2)
    data_reduced = pca.fit_transform(df_genotype_by_sample.T)
    pca_df = pd.DataFrame(data=data_reduced, columns=['PC1', 'PC2'])
    pca_df['samples'] = samples
    pca_df = pca_df.merge(pheno[['samples', 'phenotype']], on='samples', how='left')

    pca_colors = cmap(np.linspace(0, 1, len(set(pca_df.phenotype.tolist()))))
    pca_plot = os.path.join(output_dir, 'PCA.png')

    plt.figure(figsize=(10, 8))
    sns.scatterplot(
        data=pca_df.dropna(),
        x='PC1',
        y='PC2',
        hue='phenotype',
        style='phenotype',
        palette=pca_colors,
        alpha=0.7,
        s=100
    )

    plt.legend(
        title='Phenotype',
        loc='upper right',
        bbox_to_anchor=(1.4, 1),
        ncol=3
    )

    plt.title('PCA Scatter Plot')
    plt.xlabel('PCA1')
    plt.ylabel('PCA2')
    plt.grid(True)

    plt.tight_layout()

    plt.savefig(
      pca_plot,
      dpi=600
    )
    plt.close()


    # Phylogenetic tree
    dist_matrix = squareform(pdist(df_genotype_by_sample.T, metric=hamming_distance))
    linkage_matrix = hierarchy.linkage(dist_matrix, 'ward')
    ## Phylogenetic tree colors 
    hierarchy.set_link_color_palette(accent_colors)
    distance = [arr[2] for arr in linkage_matrix]
    color_threshold_q95 = np.quantile(distance, 0.95)

    tree_outfile = os.path.join(output_dir, 'tree.png')

    plt.figure(figsize=(7, 25))
    hierarchy.dendrogram(linkage_matrix, 
              labels=samples, 
              orientation='left',
              color_threshold=color_threshold_q95, 
              above_threshold_color=base_color
    )
    plt.title("Phylogenetic tree")
    plt.savefig(
      tree_outfile,
      dpi=600
    )
    plt.close()


    # GWAS
    callset = allel.read_vcf(vcf_file)
    gt = allel.GenotypeArray(callset['calldata/GT'])
    pos = callset['variants/POS']
    chrom = callset['variants/CHROM']

    keep_indices = [i for i, item in enumerate(samples) if item in pheno["samples"].tolist()]
    keep_lines = [item for i, item in enumerate(samples) if item in pheno["samples"].tolist()]

    gt_filtered = gt[:, keep_indices]
    category_dtype = pd.CategoricalDtype(categories=keep_lines, ordered=True)
    pheno['samples'] = pheno['samples'].astype(category_dtype)
    sorted_pheno = pheno.sort_values('samples')
    sorted_pheno["phenotype"] = sorted_pheno["phenotype"].astype("int")

    p_values = []
    gn = gt_filtered.to_n_alt()
    y = sorted_pheno['phenotype'].tolist()
    for i in range(gn.shape[0]):
        genotypes = gn[i]
        x = genotypes
        try:
            slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
            p_values.append(p_value)
        except Exception:
            p_values.append(np.nan)
    res_df = pd.DataFrame({'chr': chrom, 'pos': pos, 'pvalue': p_values})
    res_df['-log10(p-value)'] = -np.log10(res_df['pvalue'])
    gwas_pvals_file = os.path.join(output_dir, 'gwas_pvalues.csv')
    res_df.to_csv(gwas_pvals_file, sep="\t", index=False)

    chrom_order=natsorted(set(chrom))

    ## Adding chromosomes sizes for Manhattan plot
    chr_sizes = res_df.groupby('chr')['pos'].max().to_dict()
    chr_pos = {}
    new_pos = 0
    for Chr in chrom_order:
        chr_pos[Chr] = new_pos
        new_pos += chr_sizes[Chr]
    res_df['xpos'] = res_df.apply(lambda x: chr_pos[x.chr] + x.pos, axis=1)

    ## Manhattan plot
    mh_plot = os.path.join(output_dir, 'manhattan_plot.png')
    fig, ax = plt.subplots(figsize=(24, 6))
    sns.scatterplot(ax=ax,
                    data=res_df,
                    x='xpos',
                    y='-log10(p-value)',
                    hue='chr',
                    palette=accent_colors,
                    edgecolor='none',
                    size=100,
                    legend=False
                    )
    ticks, labels = get_xticks_and_labels(chr_pos, chr_sizes)
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    ax.set_xlabel('chromosome')
    plt.title("Manhattan plot")
    plt.savefig(
      mh_plot,
      dpi=600
    )
    plt.close()


    ## QQ-plot
    expected = -np.log10(np.linspace(0, 1, len(p_values)))
    observed = -np.log10(np.sort(p_values))

    qq_plot = os.path.join(output_dir, 'QQ_plot.png')
    plt.figure(figsize=(6, 5))
    plt.scatter(expected, observed, c=colors[0], s=10)
    plt.plot(expected, expected, color=base_color, linestyle='--')
    plt.xlabel('Expected -log10(p-value)')
    plt.ylabel('Observed -log10(p-value)')
    plt.title("QQ-plot")
    plt.savefig(
      qq_plot,
      dpi=600
    )
    plt.close()

    """
}


workflow {
    input_vcf = Channel.fromPath(params.vcf)
    input_phenotype = Channel.fromPath(params.phenotype)
    output_dir_channel = Channel.fromPath(params.output_dir)

    GWAS_analysis(input_vcf, input_phenotype, output_dir_channel)
	.view()
}
