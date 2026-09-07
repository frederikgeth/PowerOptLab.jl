"""Render standalone scientific figures from accepted study rows (matplotlib)."""
import csv, json, heapq
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
ROOT = Path(__file__).resolve().parents[2]/'studies/australian_pv'
OUT = ROOT/'figures'; OUT.mkdir(exist_ok=True)
def read(name):
    with (ROOT/'results'/name).open() as f: return list(csv.DictReader(f))
rows, devices, cases = read('customers.csv'),read('devices.csv'),read('cases.csv')
assert len(cases)==36 and all(c['accepted']=='true' for c in cases), 'Complete accepted campaign required'
assert len(rows)==36*12 and len(devices)==36*12
net=json.load((ROOT/'networks/LV3_55bus.bmopf.json').open())
labels={'unity':'Unity PF','mean_vv':'Mean VV','mean_vv_vw':'Mean VV+VW',
        'worst_vv_vw':'Worst-phase VV+VW','sequence_vv_vw':'Sequence VV+VW','sequence_droop':'Sequence + I₂ droop'}
colors=dict(zip(labels,['#555555','#d89000','#ac4875','#0072b2','#009e73','#7246a0']))
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
def save(fig,name):
    for ext in ['png','svg']:fig.savefig(OUT/f'{name}.{ext}',bbox_inches='tight')
    plt.close(fig)
def select(data,mode,level,source):
    return [r for r in data if r['mode']==mode and float(r['availability'])==level and float(r['source_pu'])==source]
def values(rs,key):return [float(r[key]) for r in rs]
# Distance follows each customer's own path; points are not joined across branches.
adj={b:[] for b in net['bus']}
for group in ['line','switch']:
    for e in net.get(group,{}).values():
        if e.get('open_switch',False):continue
        u,v=e['bus_from'],e['bus_to'];d=e.get('length',0.)
        adj[u].append((v,d));adj[v].append((u,d))
dist={'b641':0.};queue=[(0.,'b641')];parent={}
while queue:
    d,u=heapq.heappop(queue)
    if d!=dist[u]:continue
    for v,w in adj[u]:
        if d+w < dist.get(v,float('inf')):
            dist[v]=d+w;parent[v]=u;heapq.heappush(queue,(d+w,v))
customers=sorted(net['load'],key=lambda c:dist[net['load'][c]['bus']])
short={c:f'C{i+1}' for i,c in enumerate(customers)}
with (OUT/'customer_key.csv').open('w') as f:
    w=csv.writer(f, lineterminator="\n");w.writerow(['label','original_load_id','bus','path_length_m'])
    for c in customers:w.writerow([short[c],c,net['load'][c]['bus'],dist[net['load'][c]['bus']]])
fig,axs=plt.subplots(2,3,figsize=(14,7),sharex=True,sharey=True)
for ax,mode in zip(axs.flat,labels):
    rr=select(rows,mode,1.,1.03)
    for phase,color in zip('abc',['#d55e00','#0072b2','#009e73']):
        ax.scatter([dist[r['bus']] for r in rr],values(rr,f'v{phase}_pn_V'),color=color,label=phase.upper(),s=27)
    ax.axhline(253,color='black',ls='--',lw=1);ax.set_title(labels[mode]);ax.grid(alpha=.15)
    ax.set_ylim(247,257)
axs[0,0].legend(title='Phase',ncol=3,fontsize=9)
for ax in axs[1]:ax.set_xlabel('Path length from LV transformer (m)')
for ax in axs[:,0]:ax.set_ylabel('Customer phase–neutral voltage (V)')
fig.suptitle('Same feeder, loading and PV availability; only controller changes\n12 three-phase customers • 60 kW available PV • 11 kV source at 1.03 pu, 1% VUF',y=1.03)
fig.text(.5,-.02,'Upper LV supply limit: 253 V. Lower limit: 207 V (outside this zoom). Points on different branches are not connected.',ha='center')
fig.tight_layout();save(fig,'phase_voltages')
fig,axs=plt.subplots(2,3,figsize=(14,7))
for mode in labels:
    xx=[0.,.5,1.]; data=[select(rows,mode,x,1.03) for x in xx];dr=[select(devices,mode,x,1.03) for x in xx]
    ys=[ [max(values(r,'vmax_pn_V')) for r in data],
         [max(values(r,'vuf_percent')) for r in data],
         [max(values(r,'neutral_V')) for r in data],
         [sum(values(r,'poc_active_power_shortfall_W'))/1000 for r in dr],
         [-sum(values(r,'q_poc_var'))/1000 for r in dr],
         [max(values(r,'maximum_converter_current_A')) for r in dr]]
    for ax,y in zip(axs.flat,ys):ax.plot(xx,y,'o-',color=colors[mode],label=labels[mode],ms=4)
for ax,title in zip(axs.flat,['Maximum customer voltage (V)','Maximum VUF (%)','Maximum neutral displacement (V)',
                            'Available P − delivered P (kW)','Reactive power absorbed (kvar)','Maximum converter current (A)']):
    ax.set_title(title);ax.set_xlabel('PV availability fraction');ax.grid(alpha=.2)
axs[0,0].axhline(253,color='black',ls='--',lw=1)
axs[1,2].axhline(10,color='black',ls='--',lw=1)
fig.legend(*axs[0,0].get_legend_handles_labels(),loc='lower center',ncol=3,bbox_to_anchor=(.5,-.07))
fig.suptitle('Voltage benefit and inverter effort must be read together • source 1.03 pu, 1% VUF',y=1.02)
fig.tight_layout();save(fig,'tradeoffs')
# Schematic tree: y is layout only, x is route distance (not geographic position).
children={b:[] for b in dist}
for b,p in parent.items():children[p].append(b)
y={};counter=[0]
def layout(b):
    cc=sorted(children[b])
    if not cc:y[b]=counter[0];counter[0]+=1
    else:
        for c in cc:layout(c)
        y[b]=sum(y[c] for c in cc)/len(cc)
layout('b641')
fig,ax=plt.subplots(figsize=(12,6))
for b,p in parent.items():ax.plot([dist[p],dist[b]],[y[p],y[b]],color='#bbb',lw=1)
for c in customers:
    b=net['load'][c]['bus'];ax.scatter(dist[b],y[b],s=45,color='#0072b2');ax.annotate(short[c],(dist[b],y[b]),xytext=(4,4),textcoords='offset points')
ax.scatter(0,y['b641'],marker='s',s=80,color='black');ax.set_yticks([])
ax.set_xlabel('Path length from LV transformer (m)');ax.set_title('LV3: 12 customers upgraded to three-phase connections\nEach: 1.5 kW load (A/B/C = 60/25/15%) and 5 kW PV on a 6 kVA inverter')
ax.text(.01,.02,'Schematic layout; not a geographic map.\nOriginal lines, neutral, grounding and transformer retained.',transform=ax.transAxes)
save(fig,'feeder')
# Machine-readable aggregate table used by the engineering discussion.
with (ROOT/'results'/'summary.csv').open('w') as f:
    w=csv.writer(f, lineterminator="\n");w.writerow(['mode','max_V_PN','max_VUF_percent','max_neutral_V','delivered_P_kW','absorbed_Q_kvar','max_current_A'])
    for mode in labels:
        rr=select(rows,mode,1.,1.03);dd=select(devices,mode,1.,1.03)
        w.writerow([mode,max(values(rr,'vmax_pn_V')),max(values(rr,'vuf_percent')),max(values(rr,'neutral_V')),
                    sum(values(dd,'p_poc_W'))/1000,-sum(values(dd,'q_poc_var'))/1000,max(values(dd,'maximum_converter_current_A'))])
