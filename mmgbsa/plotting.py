
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from pathlib import Path
from .kinetics import calculate_ic50
import plotly.graph_objects as go
import plotly.io as pio
import plotly.express as px

def set_publication_style():
    """Set matplotlib style for publication-quality figures."""
    sns.set_style("whitegrid")
    sns.set_context("paper", font_scale=1.5)
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
    plt.rcParams['axes.linewidth'] = 1.5
    plt.rcParams['xtick.major.width'] = 1.5
    plt.rcParams['ytick.major.width'] = 1.5
    plt.rcParams['lines.linewidth'] = 2.0

def plot_violin_distributions(df, output_dir, compound_name="Ligand"):
    """
    Plot violin distributions of energy components.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing 'delta_total', 'delta_vdw', 'delta_elec', 'delta_gb', 'delta_sa' columns.
    output_dir : Path
        Directory to save the plot.
    compound_name : str
        Name of the compound.
    """
    set_publication_style()
    
    components = ['delta_total', 'delta_vdw', 'delta_elec', 'delta_gb', 'delta_sa']
    labels = ['Total', 'vdW', 'Elec', 'GB', 'SA']
    
    # Check if columns exist
    valid_components = [c for c in components if c in df.columns]
    valid_labels = [labels[i] for i, c in enumerate(components) if c in df.columns]
    
    if not valid_components:
        return None
        
    data_melted = df[valid_components].melt(var_name='Component', value_name='Energy (kcal/mol)')
    data_melted['Component'] = data_melted['Component'].replace(dict(zip(valid_components, valid_labels)))
    
    plt.figure(figsize=(10, 6))
    
    # Violin plot with inner boxplot
    sns.violinplot(x='Component', y='Energy (kcal/mol)', data=data_melted, 
                   inner='box', palette="muted", linewidth=1.5)
    
    # Add accumulation strip plot
    sns.stripplot(x='Component', y='Energy (kcal/mol)', data=data_melted, 
                  color='black', size=2, alpha=0.3, jitter=True)
    
    plt.title(f'Energy Distribution: {compound_name}')
    plt.axhline(0, color='black', linestyle='--', linewidth=1, alpha=0.5)
    
    output_file = Path(output_dir) / f"violin_distribution_{compound_name}.png"
    plt.savefig(output_file, dpi=600, bbox_inches='tight')
    plt.close()
    return output_file

def plot_convergence_enhanced(df, output_dir, compound_name="Ligand"):
    """
    Plot binding energy convergence with running mean and confidence intervals.
    """
    set_publication_style()
    
    if 'binding_energy' not in df.columns and 'delta_total' in df.columns:
        df['binding_energy'] = df['delta_total']
        
    if 'binding_energy' not in df.columns:
        return None
        
    # Calculate cumulative statistics
    cumulative_mean = df['binding_energy'].expanding().mean()
    cumulative_std = df['binding_energy'].expanding().std()
    cumulative_sem = cumulative_std / np.sqrt(np.arange(1, len(df) + 1))
    
    plt.figure(figsize=(12, 6))
    
    # Raw data (faint)
    plt.plot(df.index, df['binding_energy'], 'o', color='gray', alpha=0.1, markersize=3, label='Raw Data')
    
    # Cumulative Mean
    plt.plot(df.index, cumulative_mean, '-', color='#2980b9', linewidth=2.5, label='Cumulative Mean')
    
    # Confidence Interval (95%)
    plt.fill_between(df.index, 
                     cumulative_mean - 1.96 * cumulative_sem, 
                     cumulative_mean + 1.96 * cumulative_sem, 
                     color='#2980b9', alpha=0.2, label='95% CI')
    
    plt.xlabel('Frame')
    plt.ylabel('$\Delta G_{bind}$ (kcal/mol)')
    plt.title(f'Convergence Profile: {compound_name}')
    plt.legend(loc='best')
    
    # Add final value annotation
    final_mean = cumulative_mean.iloc[-1]
    final_sem = cumulative_sem.iloc[-1]
    plt.text(0.02, 0.05, f"Final: {final_mean:.2f} ± {final_sem:.2f} kcal/mol", 
             transform=plt.gca().transAxes, fontsize=12, 
             bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray'))
    
    output_file = Path(output_dir) / f"convergence_enhanced_{compound_name}.png"
    plt.savefig(output_file, dpi=600, bbox_inches='tight')
    plt.close()
    return output_file

def plot_res_heatmap(df, output_dir, compound_name="Ligand", top_n=20):
    """
    Plot per-residue decomposition heatmap.
    """
    set_publication_style()
    
    # Expect columns like: residue_id, total, vdw, electrostatic
    required = ['residue_id', 'total']
    if not all(c in df.columns for c in required):
        return None
        
    # Select top contributors (strongest binding = most negative total)
    df_sorted = df.sort_values('total', ascending=True).head(top_n).copy()
    
    # Prepare data for heatmap
    # Pivot isn't needed if we just want a matrix of [Residues x Components]
    # We want rows=Residues, Cols=Components
    
    cols = ['vdw', 'electrostatic', 'solvation', 'total']
    valid_cols = [c for c in cols if c in df.columns]
    
    heatmap_data = df_sorted.set_index('residue_id')[valid_cols]
    
    plt.figure(figsize=(8, len(df_sorted)*0.4 + 2))
    
    sns.heatmap(heatmap_data, annot=True, center=0, cmap="RdBu_r", fmt=".2f",
                cbar_kws={'label': 'Energy (kcal/mol)'})
    
    plt.title(f'Top {top_n} Interaction Hotspots: {compound_name}')
    plt.ylabel('Residue')
    
    output_file = Path(output_dir) / f"residue_heatmap_{compound_name}.png"
    plt.savefig(output_file, dpi=600, bbox_inches='tight')
    plt.close()
    return output_file


def plot_ic50_estimation(delta_g_mean, delta_g_sem, output_dir, compound_name="Ligand", temperature=300.0):
    """
    Visualize theoretical pIC50 with confidence intervals.
    pIC50 = -log10(IC50[M])
    """
    set_publication_style()
    
    # Calculate values
    kinetics = calculate_ic50(delta_g_mean, delta_g_sem, temperature=temperature)
    
    val = kinetics['pic50']
    low = kinetics['pic50_low']   # Low pIC50 (Weak binding)
    high = kinetics['pic50_high'] # High pIC50 (Strong binding)
    
    # Plotting
    plt.figure(figsize=(8, 4))
    
    # Horizontal Bar Plot
    # Center at val, width = high-low
    # Error bars: xerr = [[val-low], [high-val]]
    
    plt.errorbar([val], [0], xerr=[[val-low], [high-val]], fmt='o', 
                 color='black', ecolor='#8e44ad', capsize=10, elinewidth=3, markersize=10, label=f'pIC50 = {val:.2f}')
    
    plt.yticks([])
    plt.xlabel('Predicted $pIC_{50}$ (-log M)')
    plt.title(f'Predicted Potency: {compound_name}')
    
    
    # Add annotation
    plt.text(val, 0.1, f"{val:.2f}\n(95% CI: {low:.2f} - {high:.2f})", 
             ha='center', va='bottom', fontsize=14, color='#2c3e50', fontweight='bold')
             
    # Add Interpretation
    plt.text(val, -0.2, f"Theoretical $IC_{{50}}$: {kinetics['ic50']:.2e} $\mu M$", 
             ha='center', va='top', fontsize=12, color='gray')
    
    # Set limits to look nice
    margin = max(0.5, (high-low)*1.0)
    plt.xlim(low - margin, high + margin)
    plt.ylim(-0.5, 0.5)
    
    # Add grid
    plt.grid(True, axis='x', alpha=0.3)
    
    output_file = Path(output_dir) / f"ic50_prediction_{compound_name}.png"
    plt.savefig(output_file, dpi=600, bbox_inches='tight')
    plt.close()
    return output_file

def get_plotly_violin_distributions(df, compound_name="Ligand"):
    """Generate interactive Violin plot HTML."""
    components = ['delta_total', 'delta_vdw', 'delta_elec', 'delta_gb', 'delta_sa']
    labels = ['Total', 'vdW', 'Elec', 'GB', 'SA']
    
    valid_components = [c for c in components if c in df.columns]
    valid_labels = [labels[i] for i, c in enumerate(components) if c in df.columns]
    
    if not valid_components:
        return "<div>No Data for Violin Plot</div>"
        
    data_melted = df[valid_components].melt(var_name='Component', value_name='Energy')
    data_melted['Component'] = data_melted['Component'].replace(dict(zip(valid_components, valid_labels)))
    
    fig = go.Figure()
    
    for label in valid_labels:
        subset = data_melted[data_melted['Component'] == label]
        fig.add_trace(go.Violin(x=subset['Component'], y=subset['Energy'], name=label, box_visible=True, meanline_visible=True))
        
    fig.update_layout(title=f'Energy Distribution: {compound_name}', yaxis_title='Energy (kcal/mol)', template='plotly_white')
    return pio.to_html(fig, full_html=False, include_plotlyjs=False)

def get_plotly_convergence_enhanced(df, compound_name="Ligand"):
    """Generate interactive Convergence plot HTML."""
    if 'binding_energy' not in df.columns and 'delta_total' in df.columns:
        df['binding_energy'] = df['delta_total']
        
    if 'binding_energy' not in df.columns:
        return "<div>No Data for Convergence Plot</div>"
        
    cumulative_mean = df['binding_energy'].expanding().mean()
    cumulative_std = df['binding_energy'].expanding().std()
    cumulative_sem = cumulative_std / np.sqrt(np.arange(1, len(df) + 1))
    
    fig = go.Figure()
    
    # Raw Data
    fig.add_trace(go.Scatter(x=df.index, y=df['binding_energy'], mode='markers', name='Raw Data', marker=dict(color='gray', opacity=0.3, size=4)))
    
    # Cumulative Mean
    fig.add_trace(go.Scatter(x=df.index, y=cumulative_mean, mode='lines', name='Cumulative Mean', line=dict(color='#2980b9', width=3)))
    
    # Confidence Interval
    upper = cumulative_mean + 1.96 * cumulative_sem
    lower = cumulative_mean - 1.96 * cumulative_sem
    
    fig.add_trace(go.Scatter(x=df.index, y=upper, mode='lines', line=dict(width=0), showlegend=False))
    fig.add_trace(go.Scatter(x=df.index, y=lower, mode='lines', line=dict(width=0), fill='tonexty', name='95% CI', fillcolor='rgba(41, 128, 185, 0.2)'))
    
    fig.update_layout(title=f'Convergence Profile: {compound_name}', xaxis_title='Frame', yaxis_title='ΔG (kcal/mol)', template='plotly_white')
    return pio.to_html(fig, full_html=False, include_plotlyjs=False)

def get_plotly_res_heatmap(df, compound_name="Ligand", top_n=20):
    """Generate interactive Residue Heatmap HTML."""
    required = ['residue_id', 'total']
    if not all(c in df.columns for c in required):
        return "<div>No Data for Heatmap</div>"
        
    df_sorted = df.sort_values('total', ascending=True).head(top_n).copy()
    cols = ['vdw', 'electrostatic', 'solvation', 'total']
    valid_cols = [c for c in cols if c in df.columns]
    
    z = df_sorted[valid_cols].values.T
    x = df_sorted['residue_id'].tolist()
    y = valid_cols
    
    fig = go.Figure(data=go.Heatmap(z=z, x=x, y=y, colorscale='RdBu_r', zmid=0))
    fig.update_layout(title=f'Top {top_n} Interaction Hotspots: {compound_name}', xaxis_title='Residue', template='plotly_white')
    return pio.to_html(fig, full_html=False, include_plotlyjs=False)

def get_plotly_ic50_estimation(delta_g_mean, delta_g_sem, compound_name="Ligand", temperature=300.0):
    """Generate interactive pIC50 plot HTML."""
    kinetics = calculate_ic50(delta_g_mean, delta_g_sem, temperature=temperature)
    val = kinetics['pic50']
    low = kinetics['pic50_low']
    high = kinetics['pic50_high']
    
    fig = go.Figure()
    
    fig.add_trace(go.Scatter(
        x=[val], y=[0],
        error_x=dict(type='data', array=[high-val], arrayminus=[val-low], color='#8e44ad', thickness=3, width=10),
        mode='markers',
        marker=dict(color='black', size=15),
        name='pIC50'
    ))
    
    fig.update_layout(
        title=f'Predicted Potency: {compound_name}<br><sup>pIC50 = {val:.2f} (95% CI: {low:.2f} - {high:.2f})</sup>',
        xaxis_title='Predicted pIC50 (-log M)',
        yaxis=dict(showticklabels=False, range=[-0.5, 0.5]),
        template='plotly_white',
        height=300
    )
    return pio.to_html(fig, full_html=False, include_plotlyjs=False)
