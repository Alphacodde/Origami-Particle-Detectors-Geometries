import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.edgecolor'] = '#cccccc'
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['mathtext.fontset'] = 'dejavusans'
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.2), dpi=300)
geometries = ['Barrel (fine)', 'Miura-ori', 'Kresling', 'Yoshimura']
canonical_colors = ['#1f77b4', '#2ca02c', '#ff7f0e', '#d62728']
si_vals = [0.0633, 0.0646, 0.0674, 0.0688]
its3_si = 0.053
bars1 = ax1.bar(geometries, si_vals, color=canonical_colors, width=0.52, alpha=0.9, edgecolor='black', linewidth=0.5)
ax1.axhline(its3_si, color='#d95f02', linestyle='--', linewidth=1.6, label=f'ALICE ITS3 Si-only ({its3_si}%)')
for bar, val in zip(bars1, si_vals):
    ax1.text(bar.get_x() + bar.get_width() / 2, val + 0.0012, f'{val:.4f}%', ha='center', va='bottom', fontsize=8.5, fontweight='bold')
ax1.set_ylabel('Mean Radiation Length $X / X_0$ (%)', fontsize=9.5)
ax1.set_title('(a) Silicon-Only (Matched Physical Scope)', fontsize=10.5, pad=10, fontweight='bold')
ax1.set_ylim(0, 0.088)
ax1.grid(axis='y', linestyle=':', alpha=0.6)
ax1.legend(loc='upper left', frameon=True, fontsize=8.5)
ax1.tick_params(axis='x', rotation=12, labelsize=9)
stack_vals = [0.084, 0.085, 0.0868, 0.0894]
its3_full = 0.09
bars2 = ax2.bar(geometries, stack_vals, color=canonical_colors, width=0.52, alpha=0.9, edgecolor='black', linewidth=0.5)
ax2.axhline(its3_full, color='#7570b3', linestyle='--', linewidth=1.6, label=f'ITS3 full layer budget ({its3_full}%)')
for bar, val in zip(bars2, stack_vals):
    ax2.text(bar.get_x() + bar.get_width() / 2, val + 0.0015, f'{val:.4f}%', ha='center', va='bottom', fontsize=8.5, fontweight='bold')
ax2.set_ylabel('Mean Radiation Length $X / X_0$ (%)', fontsize=9.5)
ax2.set_title('(b) Si + Kapton vs. Full Engineered Layer', fontsize=10.5, pad=10, fontweight='bold')
ax2.set_ylim(0, 0.115)
ax2.grid(axis='y', linestyle=':', alpha=0.6)
ax2.legend(loc='upper left', frameon=True, fontsize=8.5)
ax2.tick_params(axis='x', rotation=12, labelsize=9)
plt.tight_layout()
out_its3 = 'Paper Template/figures/fig_its3_comparison_v2.png'
plt.savefig(out_its3, bbox_inches='tight')
print(f'Successfully generated {out_its3}')