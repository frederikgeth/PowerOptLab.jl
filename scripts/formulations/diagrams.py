"""Reproduce the documentation's scientific SVG/PDF figures.
First: julia --project scripts/formulations/figure_data.jl /tmp/pol-hinge-figure.csv
Then: python diagrams.py /tmp/pol-hinge-figure.csv
Requires numpy and matplotlib; no package runtime dependency.
"""
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
OUT=Path(__file__).resolve().parents[2]/'docs/src/assets/formulations'
OUT.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.titlesize':13,
    'axes.titleweight':'bold','axes.spines.top':False,'axes.spines.right':False,
    'axes.labelcolor':'#273548','text.color':'#273548','axes.edgecolor':'#94a3b8',
    'svg.fonttype':'none','svg.hashsalt':'PowerOptLab-formulations','savefig.facecolor':'white'})
blue,orange,green,ink='#2463a6','#cb6427','#168477','#273548'
def save(fig,name):
    fig.savefig(OUT/f'{name}.svg',bbox_inches='tight',metadata={'Date':None})
    svg=OUT/f'{name}.svg'
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
    fig.savefig(OUT/f'{name}.pdf',bbox_inches='tight',metadata={'CreationDate':None,'ModDate':None})
    fig.savefig(Path('/tmp')/f'{name}.png',dpi=120,bbox_inches='tight')
    plt.close(fig)
def panel(ax,title,xlabel,ylabel):
    ax.set(title=title,xlabel=xlabel,ylabel=ylabel)
    ax.grid(alpha=.17); ax.set_axisbelow(True)

D=np.genfromtxt(sys.argv[1],delimiter=',',names=True)
z=D['z']; delta=.05/(3/16)
fig,axes=plt.subplots(1,3,figsize=(13.8,4),layout='constrained')
for ax in axes: ax.axvspan(-delta,delta,color=green,alpha=.08)
axes[0].plot(z,D['exact'],color=ink,lw=2.7,label='Original hinge')
for field,color,label in [('soft',orange,'Softplus'),('c2',green,'Local C2'),('algebraic',blue,'Algebraic')]:
    axes[0].plot(z,D[field],color=color,lw=2,ls={'soft':'--','c2':':','algebraic':'-.'}[field],label=label)
    axes[1].plot(z,D[field+'_error'],color=color,lw=2,label=label)
    axes[2].plot(z,D[field+'_curvature'],color=color,lw=2,label=label)
panel(axes[0],'A  The function','Offset from breakpoint (V)','Hinge value (V)')
panel(axes[1],'B  Where it changes','Offset from breakpoint (V)','Approximation error (V)')
panel(axes[2],'C  Peak hinge curvature','Offset from breakpoint (V)','Second derivative (1/V)')
axes[0].legend(frameon=False,loc='upper left')
axes[1].annotate('C2 error is exactly zero\noutside the shaded patch',xy=(.48,0),xytext=(-.02,.028),
    arrowprops={'arrowstyle':'->','color':green},fontsize=10,color=green)
axes[2].text(.04,.96,'Peak = BC / error budget\nC2: 2.81, softplus: 3.47, algebraic: 5.00',
    transform=axes[2].transAxes,va='top',fontsize=9)
axes[2].set_ylim(top=6.5)
fig.suptitle('Matched maximum hinge error: 0.05 V; family constants BC determine peak curvature',fontsize=12)
save(fig,'smoothing-locality')

f=lambda v: np.clip((250-v)/10,0,1)
fig,axes=plt.subplots(1,3,figsize=(13.8,4.2),layout='constrained')
v=np.linspace(220,250,301)
axes[0].fill_between(v,0,f(v),color=green,alpha=.17)
axes[0].plot(v,f(v),color=green,lw=3)
axes[0].text(221,.12,'p ≤ 1\np ≤ (250 − V) / 10\np ≥ 0',fontsize=11)
panel(axes[0],'A  Bounded cap: convex','Voltage magnitude V (V)','Active-power fraction p')
v=np.linspace(220,270,501); roof=np.where(v<240,1,(270-v)/30)
axes[1].fill_between(v,0,f(v),color=green,alpha=.17)
axes[1].fill_between(v,f(v),roof,color=orange,alpha=.2)
axes[1].plot(v,f(v),color=green,lw=3,label='Exact cap boundary')
axes[1].plot(v,roof,'--',color=orange,lw=2,label='Hull upper boundary')
axes[1].scatter([255],[.5],color=orange,zorder=5)
axes[1].annotate('Midpoint admitted by hull\nviolates the cap',xy=(255,.5),xytext=(224,.13),
    arrowprops={'arrowstyle':'->','color':orange},fontsize=10)
panel(axes[1],'B  Floor restored: nonconvex','Voltage magnitude V (V)','Active-power fraction p')
v=np.linspace(220,250,301); chord=(250-v)/30
axes[2].fill_between(v,chord,f(v),color=orange,alpha=.2,label='Graph hull')
axes[2].plot(v,f(v),color=ink,lw=3,label='Equality graph')
axes[2].plot(v,chord,'--',color=orange,lw=2)
axes[2].scatter([245],[.08],color=green,zorder=5)
axes[2].annotate('Valid slack cap point\noutside the graph hull',xy=(245,.08),xytext=(220,.08),
    arrowprops={'arrowstyle':'->','color':green},fontsize=10)
panel(axes[2],'C  Equality ≠ upper bound','Voltage magnitude V (V)','Active-power fraction p')
axes[2].legend(frameon=False,loc='upper right')
for ax in axes: ax.set_ylim(-.04,1.12)
save(fig,'bound-geometry')

fig,ax=plt.subplots(figsize=(12.4,6.3)); ax.set(xlim=(0,12),ylim=(0,6));ax.axis('off')
def box(x,y,w,h,title,body,color):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.12,rounding_size=.1',
        facecolor=color,edgecolor='none',alpha=.12))
    ax.text(x+w/2,y+h-.28,title,ha='center',va='top',weight='bold',fontsize=12,color=color)
    ax.text(x+w/2,y+.22,body,ha='center',va='bottom',fontsize=10.5,linespacing=1.5)
def arrow(a,b): ax.annotate('',xy=b,xytext=a,arrowprops={'arrowstyle':'->','color':'#64748b','lw':1.6})
box(3.2,4.75,5.6,1.05,'ONE PHYSICAL INTENT','Volt-watt curve + sensing and power-base semantics',blue)
for x,title,body in [(0.1,'USE A: operating limit','p ≤ f(V),  220 ≤ V ≤ 250'),
                      (4.15,'USE B: operating limit','p ≤ f(V),  220 ≤ V ≤ 270'),
                      (8.2,'USE C: tracking equation','p = fδ(V),  242 ≤ V ≤ 248')]:
    box(x,2.95,3.55,1.1,title,body,blue);arrow((6,4.72),(x+1.77,4.12))
box(.1,1.05,3.55,1.25,'EXACT LINEAR ROWS','Concave cap on this domain\nNo nonlinear or integer block',green)
box(4.15,1.05,3.55,1.25,'EXPLICIT ENCODING CHOICE','Exact graph / MPCC / smoothing\nHull remains a relaxation',orange)
box(8.2,1.05,3.55,1.25,'DEPENDS ON THE FAMILY','Local C2: affine away from patches\nSoftplus: retain nonzero tails',green)
for x in (1.87,5.92,9.97): arrow((x,2.92),(x,2.36))
ax.text(6,.3,'Cache key includes the occurrence, relation, effective bounds, encoding and scales.\nEvery specialization keeps its bound assumptions as constraints.',ha='center',fontsize=11)
save(fig,'contextual-lowering')

fig,axes=plt.subplots(1,3,figsize=(13.8,4.1),layout='constrained')
v=np.linspace(230,250,301); q=-.05*(v-240); tau=.08
axes[0].fill_between(v,q-tau,q+tau,color=blue,alpha=.15)
axes[0].plot(v,q,color=ink,lw=2,label='Canonical target')
axes[0].scatter([242],[-0.06],color=orange,zorder=4)
axes[0].annotate('Optimizer may choose\nany point in this band',xy=(242,-.06),xytext=(230,-.35),
    arrowprops={'arrowstyle':'->','color':orange},fontsize=10)
panel(axes[0],'A  A set-valued model','Voltage V (V)','Reactive-power fraction q')
k=np.arange(1,9);dv=np.array([.012,.007,.003,.001,.0003,.00015,.00007,.00003]);qr=np.array([.15,.11,.06,.04,.032,.029,.027,.020])
for ax,data,tol,title,ylabel in [(axes[1],dv,1e-4,'B  Voltage settling',r'$|v_j-v_{j-1}|$ (pu)'),
                              (axes[2],qr,.025,'C  Reactive tracking',r'$|q_{desired,j}-q_{output,j}|$ (pu)')]:
    ax.semilogy(k,data,'o-',color=blue,lw=2)
    ax.axhline(tol,color=orange,ls='--',label='Illustrative tolerance')
    ax.fill_between(k,tol/5,tol,color=green,alpha=.12)
    panel(ax,title,'Control iteration j',ylabel)
    ax.legend(frameon=False,fontsize=9)
axes[1].text(.06,.08,'Voltage passes at j = 7',transform=axes[1].transAxes,fontsize=10,color=green)
axes[2].text(.06,.08,'Reactive tracking passes at j = 8',transform=axes[2].transAxes,fontsize=10,color=green)
save(fig,'tolerance-semantics')
