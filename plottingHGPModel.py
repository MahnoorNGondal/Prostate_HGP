import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from umap import UMAP
import warnings
warnings.filterwarnings('ignore')

# Set high-quality plotting style
plt.style.use('default')
sns.set_palette("husl")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 12
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['figure.facecolor'] = 'white'

print("🚀 Ultimate LPM Results Visualization Suite - FULLY FIXED")
print("="*60)

# =================================================================
# Q0: CV PERFORMANCE OVERVIEW
# =================================================================
def plot_cv_performance():
    df = pd.read_csv('q0_cv_performance.csv')
    r2_cols = ['rna', 'prot', 'ptm_combined', 'metab']
    available_cols = [col for col in r2_cols if col in df.columns]
    
    fig, axes = plt.subplots(1, 4, figsize=(20, 4))
    
    # R² scores across folds
    if len(available_cols) > 0:
        x = np.arange(1, 6)
        width = 0.2
        for i, col in enumerate(available_cols):
            values = df[col].fillna(0).clip(lower=0)
            axes[0].bar(x + i*width, values, width, alpha=0.8, 
                       label=col.replace('_', ' ').title(), edgecolor='black')
        axes[0].set_xlabel('Fold')
        axes[0].set_ylabel('R² Score')
        axes[0].set_title('5-Fold CV Performance (R²)')
        axes[0].legend()
        axes[0].set_ylim(0, max(1, max([df[col].max() for col in available_cols if col in df.columns] or [1])))
    
    # Loss trends
    axes[1].plot(df['Fold'], df['Test_Loss'], 'o-', linewidth=3, markersize=8, label='Test Loss', color='red')
    axes[1].plot(df['Fold'], df['Val_Loss'], 's--', linewidth=3, markersize=8, label='Val Loss', color='blue')
    axes[1].set_xlabel('Fold')
    axes[1].set_ylabel('Loss')
    axes[1].set_title('Reconstruction Loss')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    # Mean performance summary - FIXED
    if len(available_cols) > 0:
        summary = df[available_cols].mean().fillna(0).clip(lower=0)
        if summary.sum() > 0:
            axes[2].pie(summary.values, labels=[c.replace('_', ' ') for c in summary.index], 
                       autopct='%1.1f%%', startangle=90)
        else:
            axes[2].text(0.5, 0.5, 'No valid R²\ndata', ha='center', va='center', transform=axes[2].transAxes)
        axes[2].set_title('Mean R² Distribution')
    
    # Box plot
    df_melted = df.melt(id_vars=['Fold'], value_vars=available_cols, var_name='Modality', value_name='R2')
    df_melted['R2'] = df_melted['R2'].fillna(0).clip(lower=0)
    sns.boxplot(data=df_melted, x='Modality', y='R2', ax=axes[3])
    axes[3].set_title('R² Distribution Across Folds')
    
    plt.tight_layout()
    plt.savefig('q0_cv_performance_overview.png', bbox_inches='tight', dpi=300, facecolor='white')
    plt.close()
    print("✅ Q0: CV Performance saved")

# =================================================================
# Q1: GENE KNOCKOUT - FOXA1 (RADAR CHART FIXED)
# =================================================================
def plot_q1_gene_knockout():
    df = pd.read_csv('q1_gene_knockouts_FOXA1.csv')
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Storm intensity
    axes[0,0].bar(range(1,6), df['storm_intensity'].fillna(0), color='red', alpha=0.8, edgecolor='black')
    axes[0,0].set_xlabel('Fold')
    axes[0,0].set_ylabel('Storm Intensity')
    axes[0,0].set_title('🧬 FOXA1 Knockout: PTM Storm Intensity', fontsize=16, fontweight='bold')
    axes[0,0].grid(True, alpha=0.3)
    
    # All metrics
    metrics = ['ptm_shift', 'rna_shift', 'delta_z_norm']
    x = np.arange(1,6)
    width = 0.25
    colors = plt.cm.Set1(np.linspace(0,1,len(metrics)))
    for i, metric in enumerate(metrics):
        if metric in df.columns:
            values = df[metric].fillna(0)
            axes[0,1].bar(x + i*width, values, width, alpha=0.8, 
                         label=metric.replace('_', ' ').title(), color=colors[i])
    axes[0,1].set_xlabel('Fold')
    axes[0,1].set_ylabel('Effect Size')
    axes[0,1].set_title('Perturbation Effects')
    axes[0,1].legend()
    
    # FIXED RADAR CHART - Proper dimensions
    metrics_all = ['storm_intensity'] + [m for m in metrics if m in df.columns]
    n_metrics = len(metrics_all)
    angles = np.linspace(0, 2*np.pi, n_metrics, endpoint=False).tolist()
    means = [df[m].mean() for m in metrics_all]
    
    # Complete the circle
    angles += angles[:1]
    means += means[:1]
    
    axes[1,0].plot(angles, means, 'o-', linewidth=3, color='purple', markersize=8)
    axes[1,0].fill(angles, means, alpha=0.2, color='purple')
    axes[1,0].set_xticks(angles[:-1])
    axes[1,0].set_xticklabels([m.replace('_','\n')[:4] for m in metrics_all], fontsize=10)
    axes[1,0].set_ylim(0, max(means)*1.1)
    axes[1,0].set_title('Mean Perturbation Profile', fontsize=14, fontweight='bold')
    axes[1,0].grid(True)
    
    # Distribution
    axes[1,1].hist(df['storm_intensity'].fillna(0), bins=5, alpha=0.7, edgecolor='black', color='orange')
    axes[1,1].axvline(df['storm_intensity'].mean(), color='red', linestyle='--', linewidth=3, 
                     label=f'Mean: {df["storm_intensity"].mean():.3f}')
    axes[1,1].set_xlabel('Storm Intensity')
    axes[1,1].set_ylabel('Frequency')
    axes[1,1].set_title('Storm Intensity Distribution')
    axes[1,1].legend()
    
    plt.tight_layout()
    plt.savefig('q1_foxa1_knockout_effects.png', bbox_inches='tight', dpi=300, facecolor='white')
    plt.close()
    print("✅ Q1: FOXA1 Knockout saved")

# =================================================================
# Q2: FUSION DELETION - ERG
# =================================================================
def plot_q2_fusion_deletion():
    df = pd.read_csv('q2_fusion_deletions_ERG.csv')
    
    ptm_cols = [col for col in df.columns if any(pt in col for pt in ['pY', 'acetyl', 'methyl', 'ubiq', 'phospho'])]
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Total storm
    axes[0,0].bar(range(1,6), df['total_storm'].fillna(0), color='purple', alpha=0.8, edgecolor='black')
    axes[0,0].set_title('🔗 ERG Fusion Deletion: Total Storm', fontsize=16, fontweight='bold')
    axes[0,0].set_xlabel('Fold')
    axes[0,0].set_ylabel('Total Storm Magnitude')
    axes[0,0].grid(True, alpha=0.3)
    
    # PTM heatmap
    if ptm_cols:
        ptm_df = df[['fold'] + ptm_cols].set_index('fold').fillna(0)
        sns.heatmap(ptm_df.T, annot=True, cmap='Reds', ax=axes[0,1], 
                   cbar_kws={'label': 'Effect Size'})
        axes[0,1].set_title('PTM-Specific Effects')
    
    # Mean PTM effects
    if ptm_cols:
        ptm_means = ptm_df.mean()
        axes[1,0].bar(range(len(ptm_means)), ptm_means.values, color='skyblue', edgecolor='navy', alpha=0.8)
        axes[1,0].set_xticks(range(len(ptm_means)))
        axes[1,0].set_xticklabels(ptm_means.index, rotation=45, ha='right')
        axes[1,0].set_title('Mean PTM Disruption')
        axes[1,0].axhline(ptm_means.mean(), color='red', linestyle='--', label='Overall Mean')
        axes[1,0].legend()
    
    # Top PTMs pie
    if ptm_cols and len(ptm_cols) > 0:
        top_ptm = ptm_df.mean().nlargest(3)
        if top_ptm.sum() > 0:
            axes[1,1].pie(top_ptm.values, labels=top_ptm.index, autopct='%1.1f%%', startangle=90)
        else:
            axes[1,1].text(0.5, 0.5, 'No PTM\ndata', ha='center', va='center', transform=axes[1,1].transAxes)
    axes[1,1].set_title('Most Affected PTMs')
    
    plt.tight_layout()
    plt.savefig('q2_erg_fusion_deletion.png', bbox_inches='tight', dpi=300, facecolor='white')
    plt.close()
    print("✅ Q2: ERG Fusion Deletion saved")

# =================================================================
# Q3: MASTER REGULATORS - UMAP (PTM STORM DRIVERS)
# =================================================================
def plot_q3_master_regulators_umap():
    df = pd.read_csv('q3_master_regulators.csv')
    
    gene_means = df.groupby('gene')['saliency'].agg(['mean', 'std']).reset_index()
    gene_means['upper'] = gene_means['mean'] + gene_means['std']
    gene_means['lower'] = np.maximum(gene_means['mean'] - gene_means['std'], 0)
    
    fig = plt.figure(figsize=(20, 12))
    
    # UMAP visualization - FIXED for small datasets
    if len(gene_means) >= 2:
        try:
            reducer = UMAP(n_components=2, random_state=42, n_neighbors=min(5, len(gene_means)-1), min_dist=0.1)
            umap_coords = reducer.fit_transform(gene_means[['mean', 'std']].fillna(0).values)
            
            ax1 = plt.subplot(2, 3, 1)
            scatter = ax1.scatter(umap_coords[:,0], umap_coords[:,1], 
                                 c=gene_means['mean'].fillna(0), s=200, cmap='plasma', 
                                 alpha=0.8, edgecolors='black', linewidth=1.5)
            plt.colorbar(scatter, label='Mean Saliency')
            ax1.set_title('👑 PTM Storm Drivers (UMAP)', fontsize=16, fontweight='bold', pad=20)
            ax1.set_xlabel('UMAP 1')
            ax1.set_ylabel('UMAP 2')
        except:
            ax1 = plt.subplot(2, 3, 1)
            ax1.text(0.5, 0.5, 'UMAP\nNot Available', ha='center', va='center', transform=ax1.transAxes)
            ax1.set_title('👑 PTM Storm Drivers (UMAP)')
    
    # Top 10 genes
    top_genes = gene_means.nlargest(10, 'mean')
    ax2 = plt.subplot(2, 3, 2)
    colors = plt.cm.plasma(np.linspace(0,1,len(top_genes)))
    bars = ax2.barh(range(len(top_genes)), top_genes['mean'], 
                   xerr=[top_genes['mean']-top_genes['lower'], top_genes['upper']-top_genes['mean']],
                   capsize=5, color=colors, alpha=0.8, edgecolor='black')
    ax2.set_yticks(range(len(top_genes)))
    ax2.set_yticklabels(top_genes['gene'])
    ax2.set_xlabel('Saliency Score')
    ax2.set_title('Top 10 PTM Storm Drivers')
    
    # All genes ranked
    ax3 = plt.subplot(2, 3, 3)
    ax3.bar(range(len(gene_means)), gene_means['mean'].fillna(0), 
           yerr=gene_means['std'].fillna(0), capsize=3, alpha=0.7, 
           color='steelblue', edgecolor='navy')
    ax3.axhline(y=gene_means['mean'].mean(), color='red', linestyle='--', label='Overall Mean')
    ax3.set_ylabel('Saliency Score')
    ax3.set_title('All Genomic Drivers Ranked')
    ax3.legend()
    
    # Cumulative distribution
    ax4 = plt.subplot(2, 3, 4)
    sorted_means = np.sort(gene_means['mean'].fillna(0))[::-1]
    ax4.plot(np.cumsum(sorted_means)/np.sum(sorted_means), 'o-', linewidth=3, markersize=6)
    ax4.axhline(0.8, color='red', linestyle='--', label='80% Storm Power')
    ax4.set_xlabel('Top N Genes')
    ax4.set_ylabel('Cumulative Storm Contribution')
    ax4.set_title('Pareto Analysis: 80/20 Rule')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    # Consistency heatmap
    ax5 = plt.subplot(2, 3, 5)
    pivot = df.pivot(index='gene', columns='fold', values='saliency').fillna(0)
    top15_genes = gene_means.nlargest(15, 'mean')['gene'].tolist()
    if len(top15_genes) > 0 and len(pivot) > 0:
        try:
            sns.heatmap(pivot.loc[top15_genes], annot=True, cmap='YlOrRd', ax=ax5)
        except:
            sns.heatmap(pivot.head(15), annot=True, cmap='YlOrRd', ax=ax5)
    ax5.set_title('Consistency Across Folds (Top 15)')
    
    plt.tight_layout()
    plt.savefig('q3_ptm_storm_drivers_umap.png', bbox_inches='tight', dpi=300, facecolor='white')
    plt.close()
    print("✅ Q3: PTM Storm Drivers UMAP saved")

# =================================================================
# Q4-Q6: SIMPLIFIED & ROBUST
# =================================================================
def plot_q4_resistance_drivers():
    df = pd.read_csv('q4_resistance_drivers.csv')
    gene_means = df.groupby('gene')['resistance_score'].agg(['mean', 'std']).sort_values('mean', ascending=False).head(10)
    
    fig, ax = plt.subplots(figsize=(12, 8))
    bars = ax.barh(range(len(gene_means)), gene_means['mean'], 
                   xerr=gene_means['std'], capsize=5, color='darkred', alpha=0.8, edgecolor='black')
    ax.set_yticks(range(len(gene_means)))
    ax.set_yticklabels(gene_means.index)
    ax.set_xlabel('Resistance Score')
    ax.set_title('🛡️ Top 10 Resistance Drivers', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig('q4_resistance_drivers.png', bbox_inches='tight', dpi=300, facecolor='white')
    plt.close()
    print("✅ Q4: Resistance Drivers saved")

def plot_q5_ptm_dependencies():
    df = pd.read_csv('q5_ptm_dependencies.csv')
    ptm_stats = df.groupby('ptm')['dependency_strength'].mean().sort_values(ascending=False)
    
    fig, ax = plt.subplots(figsize=(12, 8))
    bars = ax.bar(range(len(ptm_stats)), ptm_stats.values, color='purple', alpha=0.8, edgecolor='black')
    ax.set_xticks(range(len(ptm_stats)))
    ax.set_xticklabels(ptm_stats.index, rotation=45, ha='right')
    ax.set_ylabel('Dependency Strength')
    ax.set_title('🔗 PTM Dependency Strength', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig('q5_ptm_dependencies.png', bbox_inches='tight', dpi=300, facecolor='white')
    plt.close()
    print("✅ Q5: PTM Dependencies saved")

def plot_q6_storm_leaders():
    df = pd.read_csv('q6_storm_leaders.csv')
    ptm_stats = df.groupby('ptm')['metabolome_r2'].mean().sort_values(ascending=False)
    
    fig, ax = plt.subplots(figsize=(12, 8))
    bars = ax.barh(range(len(ptm_stats)), ptm_stats.values, color='gold', alpha=0.9, edgecolor='darkgoldenrod')
    ax.set_yticks(range(len(ptm_stats)))
    ax.set_yticklabels(ptm_stats.index)
    ax.set_xlabel('Metabolome R²')
    ax.set_title('⭐ PTM Storm Leaders', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig('q6_storm_leaders.png', bbox_inches='tight', dpi=300, facecolor='white')
    plt.close()
    print("✅ Q6: Storm Leaders saved")

# =================================================================
# EXECUTE ALL
# =================================================================
if __name__ == "__main__":
    print("\n🎨 Generating Slide-Ready Visualizations...\n")
    
    plot_cv_performance()
    plot_q1_gene_knockout()
    plot_q2_fusion_deletion()
    plot_q3_master_regulators_umap()
    plot_q4_resistance_drivers()
    plot_q5_ptm_dependencies()
    plot_q6_storm_leaders()
    
    print("\n" + "="*60)
    print("🏆 ALL 7 SLIDE-READY PLOTS GENERATED SUCCESSFULLY!")
    print("📁 Files saved:")
    print("- q0_cv_performance_overview.png")
    print("- q1_foxa1_knockout_effects.png")
    print("- q2_erg_fusion_deletion.png") 
    print("- q3_ptm_storm_drivers_umap.png")
    print("- q4_resistance_drivers.png")
    print("- q5_ptm_dependencies.png")
    print("- q6_storm_leaders.png")
    print("\n🎯 READY FOR PRESENTATION! 300 DPI publication quality.")
