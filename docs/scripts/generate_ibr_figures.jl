#!/usr/bin/env julia

# Equation-generated figures for the IBR technical chapter. This script uses
# only Julia standard libraries so documentation figures remain reproducible in
# the docs environment without adding a plotting dependency.

using Printf

const OUT = normpath(joinpath(@__DIR__, "..", "src", "assets", "ibr"))

points(xs, ys) = join((@sprintf("%.2f,%.2f", x, y) for (x, y) in zip(xs, ys)), " ")

function waveform_figure(path)
    a = cis(2pi/3)
    vph = 416.0 / sqrt(3)
    iref = 20.0 / sqrt(2)
    V = vph .* ComplexF64[1, a^-1, a^-2]
    I = iref .* ComplexF64[1, a, a^2] # Deakin et al. Table V, case 3d
    S2 = sum(V .* I)
    vdc, cdc, f = 700.0, 50e-3, 50.0
    denom = 2 * (2pi*f) * cdc * vdc
    D = im * S2 / denom

    theta = range(0, 2pi; length=721)
    time_ms = collect(theta) ./ (2pi*f) .* 1e3
    p2 = [real(S2 * cis(2th)) / 1e3 for th in theta]
    dv = [real(D * cis(2th)) for th in theta]

    x0, x1 = 95.0, 955.0
    mapx(t) = x0 + (x1-x0) * t / maximum(time_ms)
    mapy1(y) = 180.0 - 100.0 * y / (abs(S2)/1e3)
    mapy2(y) = 435.0 - 95.0 * y / abs(D)
    xs = mapx.(time_ms)
    ppts = points(xs, mapy1.(p2))
    dpts = points(xs, mapy2.(dv))
    gridx = join(("<line x1=\"$(mapx(t))\" y1=\"70\" x2=\"$(mapx(t))\" y2=\"535\"/>"
                  for t in 0:5:20), "\n")

    svg = """<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="620" viewBox="0 0 1000 620">
<rect width="1000" height="620" fill="#ffffff"/>
<style>
 text { font-family: Inter, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; fill:#172033 }
 .grid { stroke:#dfe5ed; stroke-width:1 }
 .axis { stroke:#526071; stroke-width:1.4 }
 .small { font-size:14px; fill:#526071 }
 .label { font-size:17px; font-weight:600 }
</style>
<text x="50" y="30" font-size="23" font-weight="700">Published negative-sequence case: AC 2ω power and DC-link response</text>
<text class="small" x="50" y="55">Table V case 3d · 416 V L-L RMS · 20 A peak · Vdc=700 V · Cdc=50 mF · f=50 Hz</text>
<g class="grid">$gridx
 <line x1="$x0" y1="80" x2="$x1" y2="80"/><line x1="$x0" y1="180" x2="$x1" y2="180"/><line x1="$x0" y1="280" x2="$x1" y2="280"/>
 <line x1="$x0" y1="340" x2="$x1" y2="340"/><line x1="$x0" y1="435" x2="$x1" y2="435"/><line x1="$x0" y1="530" x2="$x1" y2="530"/>
</g>
<line class="axis" x1="$x0" y1="180" x2="$x1" y2="180"/>
<line class="axis" x1="$x0" y1="435" x2="$x1" y2="435"/>
<polyline points="$ppts" fill="none" stroke="#1769aa" stroke-width="3"/>
<polyline points="$dpts" fill="none" stroke="#d65a31" stroke-width="3"/>
<text class="label" x="25" y="101">p₂ω</text><text class="small" x="25" y="120">(kW)</text>
<text class="small" x="35" y="84">+$(@sprintf("%.2f", abs(S2)/1e3))</text>
<text class="small" x="52" y="185">0</text>
<text class="small" x="33" y="284">−$(@sprintf("%.2f", abs(S2)/1e3))</text>
<text class="label" x="25" y="365">Δvdc</text><text class="small" x="25" y="384">(V peak)</text>
<text class="small" x="35" y="344">+$(@sprintf("%.3f", abs(D)))</text>
<text class="small" x="52" y="440">0</text>
<text class="small" x="31" y="534">−$(@sprintf("%.3f", abs(D)))</text>
<g class="small">$(join(("<text x=\"$(mapx(t)-6)\" y=\"561\">$t</text>" for t in 0:5:20), ""))</g>
<text class="label" x="472" y="589">time (ms)</text>
<text class="small" x="95" y="613">|S̃|=$(@sprintf("%.3f", abs(S2)/1e3)) kW; |D|=$(@sprintf("%.3f", abs(D))) V peak. D=jS̃/(2ωCdcVdc), so capacitor energy balance introduces a 90° shift.</text>
</svg>"""
    write(path, svg)
end

function utilisation_figure(path)
    m = 0.96
    U = collect(range(0.0, 280.0; length=281))
    v_three = sqrt(6) .* U ./ m
    v_split = 2sqrt(2) .* U ./ m
    x0, x1, y0, y1 = 92.0, 950.0, 510.0, 72.0
    mapx(u) = x0 + (x1-x0)*u/280
    mapy(v) = y0 - (y0-y1)*v/900
    pthree = points(mapx.(U), mapy.(v_three))
    psplit = points(mapx.(U), mapy.(v_split))
    uref = 230.0
    v3 = sqrt(6)*uref/m
    vs = 2sqrt(2)*uref/m
    xref = mapx(uref)
    hgrid = join(("<line x1=\"$x0\" y1=\"$(mapy(v))\" x2=\"$x1\" y2=\"$(mapy(v))\"/>"
                  for v in 0:150:900), "\n")
    vgrid = join(("<line x1=\"$(mapx(u))\" y1=\"$y1\" x2=\"$(mapx(u))\" y2=\"$y0\"/>"
                  for u in 0:40:280), "\n")

    svg = """<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="610" viewBox="0 0 1000 610">
<rect width="1000" height="610" fill="#ffffff"/>
<style>
 text { font-family: Inter, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; fill:#172033 }
 .grid { stroke:#dfe5ed; stroke-width:1 }
 .axis { stroke:#526071; stroke-width:1.5 }
 .small { font-size:14px; fill:#526071 }
 .label { font-size:17px; font-weight:600 }
</style>
<text x="50" y="38" font-size="23" font-weight="700">Balanced fundamental voltage utilisation of the ideal switching hull</text>
<g class="grid">$hgrid$vgrid</g>
<line class="axis" x1="$x0" y1="$y0" x2="$x1" y2="$y0"/>
<line class="axis" x1="$x0" y1="$y0" x2="$x0" y2="$y1"/>
<polyline points="$pthree" fill="none" stroke="#1769aa" stroke-width="4"/>
<polyline points="$psplit" fill="none" stroke="#d65a31" stroke-width="4"/>
<line x1="$xref" y1="$y0" x2="$xref" y2="$(mapy(vs))" stroke="#64748b" stroke-width="1.5" stroke-dasharray="7 6"/>
<circle cx="$xref" cy="$(mapy(v3))" r="6" fill="#1769aa"/>
<circle cx="$xref" cy="$(mapy(vs))" r="6" fill="#d65a31"/>
<text class="label" x="565" y="270" fill="#1769aa">3-leg / 4-leg pairwise hull</text>
<text class="small" x="565" y="291">Vdc = √6 U/m</text>
<text class="label" x="565" y="155" fill="#d65a31">split DC link</text>
<text class="small" x="565" y="176">Vdc = 2√2 U/m</text>
<rect x="690" y="394" width="250" height="96" rx="8" fill="#f4f7fa" stroke="#ccd6e2"/>
<text class="small" x="706" y="419">At U = 230 V RMS and m = 0.96:</text>
<text class="small" x="706" y="444">3-/4-leg: $(@sprintf("%.1f", v3)) V</text>
<text class="small" x="706" y="466">split link: $(@sprintf("%.1f", vs)) V (+15.5%)</text>
$(join(("<text class=\"small\" x=\"$(mapx(u)-8)\" y=\"535\">$u</text>" for u in 0:40:280), ""))
$(join(("<text class=\"small\" x=\"48\" y=\"$(mapy(v)+5)\">$v</text>" for v in 0:150:900), ""))
<text class="label" x="410" y="574">phase voltage U (V RMS)</text>
<text class="label" transform="translate(20 390) rotate(-90)">minimum total DC voltage (V)</text>
<text class="small" x="93" y="598">Balanced, zero-ripple, ideal two-level boundary. The 4-leg per-phase rail is non-binding for a balanced set.</text>
</svg>"""
    write(path, svg)
end

function lcl_circuit_figure(path)
    svg = """<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="430" viewBox="0 0 1000 430">
<rect width="1000" height="430" fill="#ffffff"/>
<style>
 text { font-family: Inter, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; fill:#172033 }
 .title { font-size:23px; font-weight:700 }
 .label { font-size:17px; font-weight:600 }
 .small { font-size:14px; fill:#526071 }
 .wire { stroke:#334155; stroke-width:3; fill:none }
 .arrow { stroke:#1769aa; stroke-width:2.5; fill:none; marker-end:url(#arrow) }
 .box { fill:#f4f7fa; stroke:#8ca0b5; stroke-width:2 }
</style>
<defs><marker id="arrow" markerWidth="8" markerHeight="8" refX="7" refY="3" orient="auto"><path d="M0,0 L0,6 L8,3 z" fill="#1769aa"/></marker></defs>
<text class="title" x="42" y="38">Explicit fundamental-frequency LCL circuit and reference directions</text>

<rect class="box" x="48" y="128" width="145" height="112" rx="8"/>
<text class="label" x="78" y="170">converter</text>
<text class="small" x="67" y="196">internal voltage Ux</text>
<text class="small" x="67" y="218">switching hull + loss</text>

<line class="wire" x1="193" y1="184" x2="273" y2="184"/>
<rect class="box" x="273" y="140" width="148" height="88" rx="8"/>
<text class="label" x="303" y="174">converter arm</text>
<text class="small" x="300" y="202">primitive Zc = Rc+jXc</text>
<line class="wire" x1="421" y1="184" x2="512" y2="184"/>
<circle cx="512" cy="184" r="7" fill="#172033"/>
<text class="label" x="465" y="112">midpoint Mx</text>
<text class="small" x="455" y="132">reported Vfilter</text>

<line class="wire" x1="519" y1="184" x2="602" y2="184"/>
<rect class="box" x="602" y="140" width="145" height="88" rx="8"/>
<text class="label" x="635" y="174">grid arm</text>
<text class="small" x="623" y="202">primitive Zg = Rg+jXg</text>
<line class="wire" x1="747" y1="184" x2="824" y2="184"/>
<circle cx="824" cy="184" r="7" fill="#172033"/>
<rect class="box" x="824" y="128" width="128" height="112" rx="8"/>
<text class="label" x="857" y="170">POC</text>
<text class="small" x="844" y="196">network Vp,xn</text>
<text class="small" x="844" y="218">grid current Ig</text>

<line class="arrow" x1="215" y1="78" x2="391" y2="78"/>
<text class="small" x="270" y="69">converter current Ic</text>
<line class="arrow" x1="539" y1="78" x2="716" y2="78"/>
<text class="small" x="595" y="69">grid current Ig</text>

<line class="wire" x1="512" y1="191" x2="512" y2="252"/>
<rect class="box" x="466" y="252" width="92" height="45" rx="6"/>
<text class="small" x="483" y="280">Rd damping</text>
<line class="wire" x1="512" y1="297" x2="512" y2="322"/>
<line class="wire" x1="484" y1="322" x2="540" y2="322"/>
<line class="wire" x1="484" y1="337" x2="540" y2="337"/>
<line class="wire" x1="512" y1="337" x2="512" y2="369"/>
<line class="wire" x1="426" y1="369" x2="598" y2="369"/>
<text class="label" x="558" y="329">Cmid</text>
<text class="small" x="609" y="374">midpoint neutral / reference</text>
<line class="arrow" x1="444" y1="238" x2="444" y2="315"/>
<text class="small" x="351" y="278">Ish=(G+jB)M</text>

<rect x="57" y="312" width="315" height="76" rx="8" fill="#eef6fb" stroke="#a8c9df"/>
<text class="small" x="75" y="337">KVL: U-M = Zc Ic; M-Vp = Zg Ig</text>
<text class="small" x="75" y="361">KCL: Ic = Ig + Ish</text>
<text class="small" x="75" y="383">Ysh = 1/(Rd + 1/(j omega Cmid))</text>

<text class="small" x="42" y="416">Each drawn phase has a phase-to-neutral equation; both series arms may instead use full conductor-domain primitive matrices.</text>
</svg>"""
    write(path, svg)
end

function split_link_figure(path)
    svg = """<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="610" viewBox="0 0 1000 610">
<rect width="1000" height="610" fill="#ffffff"/>
<style>
 text { font-family: Inter, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; fill:#172033 }
 .title { font-size:23px; font-weight:700 }
 .label { font-size:17px; font-weight:600 }
 .small { font-size:14px; fill:#526071 }
 .wire { stroke:#334155; stroke-width:3; fill:none }
 .component { fill:#f4f7fa; stroke:#8ca0b5; stroke-width:2 }
 .neutral { stroke:#1769aa; stroke-width:3; fill:none; marker-end:url(#blue-arrow) }
 .common { stroke:#d65a31; stroke-width:3; fill:none; marker-end:url(#orange-arrow) }
 .balance { stroke:#33865a; stroke-width:3; fill:none; marker-start:url(#green-back); marker-end:url(#green-arrow) }
 .formula { fill:#f7f9fc; stroke:#ccd6e2; stroke-width:1.5 }
</style>
<defs>
 <marker id="blue-arrow" markerWidth="8" markerHeight="8" refX="7" refY="3" orient="auto"><path d="M0,0 L0,6 L8,3 z" fill="#1769aa"/></marker>
 <marker id="orange-arrow" markerWidth="8" markerHeight="8" refX="7" refY="3" orient="auto"><path d="M0,0 L0,6 L8,3 z" fill="#d65a31"/></marker>
 <marker id="green-arrow" markerWidth="8" markerHeight="8" refX="7" refY="3" orient="auto"><path d="M0,0 L0,6 L8,3 z" fill="#33865a"/></marker>
 <marker id="green-back" markerWidth="8" markerHeight="8" refX="1" refY="3" orient="auto"><path d="M8,0 L8,6 L0,3 z" fill="#33865a"/></marker>
</defs>
<text class="title" x="42" y="38">Asymmetric split DC link: charge, current sharing, and thermal stress</text>
<text class="small" x="42" y="62">Three-leg, four-wire bridge · physical RMS and thermally equivalent currents are reported separately</text>

<line class="wire" x1="92" y1="102" x2="682" y2="102"/>
<line class="wire" x1="92" y1="478" x2="682" y2="478"/>
<text class="label" x="50" y="108">+Vdc/2</text>
<text class="label" x="50" y="484">−Vdc/2</text>

<line class="wire" x1="310" y1="102" x2="310" y2="137"/>
<rect class="component" x="270" y="137" width="80" height="38" rx="4"/>
<text class="small" x="279" y="162">ESR upper</text>
<line class="wire" x1="310" y1="175" x2="310" y2="213"/>
<line class="wire" x1="274" y1="213" x2="346" y2="213"/>
<line class="wire" x1="274" y1="229" x2="346" y2="229"/>
<line class="wire" x1="310" y1="229" x2="310" y2="290"/>
<text class="label" x="360" y="227">Cu</text>

<circle cx="310" cy="290" r="7" fill="#172033"/>
<text class="label" x="270" y="278">midpoint M</text>

<line class="wire" x1="310" y1="290" x2="310" y2="351"/>
<line class="wire" x1="274" y1="351" x2="346" y2="351"/>
<line class="wire" x1="274" y1="367" x2="346" y2="367"/>
<line class="wire" x1="310" y1="367" x2="310" y2="405"/>
<rect class="component" x="270" y="405" width="80" height="38" rx="4"/>
<text class="small" x="279" y="430">ESR lower</text>
<line class="wire" x1="310" y1="443" x2="310" y2="478"/>
<text class="label" x="360" y="365">Cl</text>

<rect class="component" x="565" y="151" width="220" height="194" rx="9"/>
<text class="label" x="606" y="180">three phase legs</text>
<text class="small" x="596" y="205">each leg spans the total bus</text>
<line class="wire" x1="682" y1="102" x2="682" y2="151"/>
<line class="wire" x1="682" y1="345" x2="682" y2="478"/>
<line class="wire" x1="785" y1="220" x2="855" y2="220"/>
<line class="wire" x1="785" y1="260" x2="855" y2="260"/>
<line class="wire" x1="785" y1="300" x2="855" y2="300"/>
<text class="label" x="866" y="226">a</text><text class="label" x="866" y="266">b</text><text class="label" x="866" y="306">c</text>
<line class="wire" x1="317" y1="290" x2="448" y2="290"/>
<line class="wire" x1="448" y1="290" x2="448" y2="380"/>
<line x1="663" y1="380" x2="701" y2="380" stroke="#ffffff" stroke-width="11"/>
<line class="wire" x1="448" y1="380" x2="855" y2="380"/>
<text class="label" x="866" y="386">n</text>

<line class="neutral" x1="538" y1="270" x2="351" y2="270"/>
<text class="small" x="398" y="259" fill="#1769aa">neutral In</text>
<text class="small" x="111" y="194" fill="#1769aa">upper share αu In</text>
<line class="neutral" x1="248" y1="190" x2="248" y2="249"/>
<text class="small" x="111" y="397" fill="#1769aa">lower share αl In</text>
<line class="neutral" x1="248" y1="390" x2="248" y2="331"/>

<line class="common" x1="393" y1="142" x2="393" y2="235"/>
<line class="common" x1="393" y1="345" x2="393" y2="438"/>
<text class="small" x="408" y="157" fill="#d65a31">common I2ω + Isw</text>
<text class="small" x="408" y="424" fill="#d65a31">common I2ω + Isw</text>

<line class="balance" x1="520" y1="190" x2="520" y2="360"/>
<text class="small" x="532" y="358" fill="#33865a">bounded charge transfer qb</text>

<rect class="formula" x="48" y="516" width="904" height="70" rx="8"/>
<text class="small" x="65" y="540">Ceq = Cu Cl/(Cu+Cl)     αu = Cu/(Cu+Cl), αl = Cl/(Cu+Cl)     |N| = |In|/[ω(Cu+Cl)]</text>
<text class="small" x="65" y="565">Nmean = Vdc(Cu−Cl)/[2(Cu+Cl)] + 2qb/(Cu+Cl)     Ih² = αh²|In|² + I2ω² + Isw²</text>
</svg>"""
    write(path, svg)
end

function pwm_ripple_figure(path)
    svg = """<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="640" viewBox="0 0 1000 640">
<rect width="1000" height="640" fill="#ffffff"/>
<style>
 .title { font: 700 24px sans-serif; fill:#172033 }
 .subtitle { font: 14px sans-serif; fill:#526174 }
 .head { font: 700 16px sans-serif; fill:#172033 }
 .body { font: 13px sans-serif; fill:#334155 }
 .formula { font: 14px serif; fill:#172033 }
 .box { fill:#f7f9fc; stroke:#8ca0b8; stroke-width:1.7 }
 .audit { fill:#edf7ff; stroke:#1769aa; stroke-width:2 }
 .constraint { fill:#fff6eb; stroke:#d65a31; stroke-width:2 }
 .result { fill:#edf8f1; stroke:#33865a; stroke-width:2 }
 .arrow { stroke:#1769aa; stroke-width:2.3; fill:none; marker-end:url(#arrow) }
 .feedback { stroke:#33865a; stroke-width:2.6; fill:none; marker-end:url(#feedback) }
 .dash { stroke:#94a3b8; stroke-width:1.5; stroke-dasharray:6 5; fill:none }
</style>
<defs>
 <marker id="arrow" markerWidth="9" markerHeight="9" refX="8" refY="3" orient="auto"><path d="M0,0 L0,6 L9,3 z" fill="#1769aa"/></marker>
 <marker id="feedback" markerWidth="9" markerHeight="9" refX="8" refY="3" orient="auto"><path d="M0,0 L0,6 L9,3 z" fill="#33865a"/></marker>
</defs>
<text class="title" x="42" y="40">Carrier-level PWM ripple audit and conservative closure</text>
<text class="subtitle" x="42" y="65">Fundamental phasors remain optimisation variables; ideal shared-carrier states provide an operating-point oracle</text>

<rect class="box" x="45" y="105" width="230" height="174" rx="10"/>
<text class="head" x="67" y="136">1 · Smooth phasor NLP</text>
<text class="body" x="67" y="166">solves U, I, Vdc and 2ω ripple</text>
<text class="body" x="67" y="190">subject to topology, filter,</text>
<text class="body" x="67" y="214">sequence and capacitor bounds</text>
<text class="formula" x="67" y="250">Ialloc = κ √Σ |Ileg|²</text>

<path class="arrow" d="M275 192 L345 192"/>
<text class="body" x="288" y="178">U, I, Vdc</text>

<rect class="audit" x="345" y="92" width="305" height="285" rx="10"/>
<text class="head" x="368" y="124">2 · Ideal shared-carrier audit</text>
<text class="body" x="368" y="154">SPWM: γ = 0</text>
<text class="body" x="368" y="177">centered PWM: γ = −(max u + min u)/2</text>
<text class="formula" x="368" y="211">dℓ = 1/2 + (uℓ + γ)/vdc</text>
<rect x="375" y="233" width="245" height="48" rx="5" fill="#ffffff" stroke="#77a9d1"/>
<polyline points="390,265 420,241 450,265 480,241 510,265 540,241 570,265 600,241" fill="none" stroke="#1769aa" stroke-width="2"/>
<text class="body" x="389" y="303">one triangular carrier → correlated sℓ</text>
<text class="formula" x="389" y="330">î = Σ (sℓ − dℓ) iℓ</text>
<text class="body" x="389" y="355">integrate î² and capacitor charge</text>

<path class="arrow" d="M650 192 L720 192"/>
<text class="body" x="664" y="178">audit</text>

<rect class="result" x="720" y="105" width="235" height="220" rx="10"/>
<text class="head" x="742" y="136">3 · Physical diagnostics</text>
<text class="body" x="742" y="168">predicted switching RMS Isw</text>
<text class="body" x="742" y="193">voltage ripple Vrms and Vpp</text>
<text class="body" x="742" y="218">minimum duty headroom</text>
<line x1="742" y1="240" x2="932" y2="240" stroke="#9bc7aa"/>
<text class="formula" x="742" y="268">ΔI = Ialloc − Isw</text>
<text class="body" x="742" y="295">publish only when ΔI ≥ 0</text>

<rect class="constraint" x="95" y="445" width="810" height="108" rx="10"/>
<text class="head" x="120" y="477">Capacitor constraint retained inside every smooth solve</text>
<text class="formula" x="120" y="510">wn In² + w2 I2ω² + wsw (i_sw² + κ² Σ |Ileg|²) ≤ Irated²</text>
<text class="body" x="120" y="535">manual i_sw is an independent residual · upper/lower split banks are constrained separately</text>

<path class="feedback" d="M837 325 L837 405 Q837 425 817 425 L183 425 Q160 425 160 445"/>
<text class="body" x="561" y="416" fill="#33865a">update κ from the audited operating point</text>
<path class="feedback" d="M95 500 L60 500 Q35 500 35 475 L35 230 Q35 205 45 205"/>
<text class="body" x="44" y="592">Repeat until the allocated reserve covers the prediction within tolerance.</text>
<path class="dash" d="M345 590 L650 590"/>
<text class="subtitle" x="674" y="595">not switching-network EMT</text>
</svg>"""
    write(path, svg)
end

function ac_ripple_network_figure(path)
    svg = """<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="620" viewBox="0 0 1000 620">
<rect width="1000" height="620" fill="#ffffff"/>
<style>
 .title { font:700 24px sans-serif; fill:#172033 }
 .sub { font:14px sans-serif; fill:#526174 }
 .head { font:700 16px sans-serif; fill:#172033 }
 .body { font:13px sans-serif; fill:#334155 }
 .formula { font:14px serif; fill:#172033 }
 .wire { stroke:#172033; stroke-width:2.4; fill:none }
 .harm { stroke:#1769aa; stroke-width:2.5; fill:none; marker-end:url(#blue) }
 .current { stroke:#d65a31; stroke-width:2.4; fill:none; marker-end:url(#orange) }
 .box { fill:#f7f9fc; stroke:#8ca0b8; stroke-width:1.6 }
 .bluebox { fill:#edf7ff; stroke:#1769aa; stroke-width:1.8 }
 .greenbox { fill:#edf8f1; stroke:#33865a; stroke-width:1.8 }
</style>
<defs>
 <marker id="blue" markerWidth="9" markerHeight="9" refX="8" refY="3" orient="auto"><path d="M0,0 L0,6 L9,3 z" fill="#1769aa"/></marker>
 <marker id="orange" markerWidth="9" markerHeight="9" refX="8" refY="3" orient="auto"><path d="M0,0 L0,6 L9,3 z" fill="#d65a31"/></marker>
</defs>
<text class="title" x="42" y="40">AC switching ripple: carrier harmonics through the physical filter</text>
<text class="sub" x="42" y="65">Every carrier harmonic sees the conductor primitive, neutral return, LCL midpoint branch, and grid-side arm</text>

<rect class="bluebox" x="45" y="115" width="180" height="155" rx="10"/>
<text class="head" x="65" y="145">shared-carrier bridge</text>
<text class="body" x="65" y="176">3-leg: sum-zero projection</text>
<text class="body" x="65" y="199">split: physical midpoint</text>
<text class="body" x="65" y="222">4-leg: phase − neutral pole</text>
<text class="formula" x="65" y="250">Êh = F{Vdc(s−d)}</text>

<path class="harm" d="M225 192 L315 192"/>
<text class="body" x="246" y="178">Êh</text>

<rect class="box" x="315" y="135" width="160" height="115" rx="8"/>
<text class="head" x="340" y="165">converter arm</text>
<text class="formula" x="337" y="197">Zc,h = Rc + jΩhLc</text>
<text class="body" x="337" y="225">primitive + neutral/mutual</text>
<path class="wire" d="M475 192 L570 192"/>
<text class="body" x="492" y="178">midpoint Vm,h</text>

<circle cx="570" cy="192" r="6" fill="#172033"/>
<path class="wire" d="M570 192 L570 340"/>
<rect class="box" x="500" y="340" width="140" height="92" rx="8"/>
<text class="head" x="523" y="369">damped shunt</text>
<text class="formula" x="520" y="400">Rd + 1/(jΩhCf)</text>
<path class="wire" d="M570 432 L570 468"/>
<line x1="535" y1="468" x2="605" y2="468" stroke="#172033" stroke-width="2.4"/>
<line x1="545" y1="478" x2="595" y2="478" stroke="#172033" stroke-width="2.4"/>
<line x1="555" y1="488" x2="585" y2="488" stroke="#172033" stroke-width="2.4"/>

<path class="wire" d="M576 192 L665 192"/>
<rect class="box" x="665" y="135" width="160" height="115" rx="8"/>
<text class="head" x="701" y="165">grid arm</text>
<text class="formula" x="687" y="197">Zg,h = Rg + jΩhLg</text>
<text class="body" x="687" y="225">primitive + neutral/mutual</text>
<path class="wire" d="M825 192 L935 192"/>
<rect x="935" y="154" width="8" height="76" fill="#172033"/>
<text class="body" x="866" y="178">stiff POC</text>

<path class="current" d="M445 285 L335 285"/>
<text class="body" x="348" y="307" fill="#d65a31">converter Ic,h</text>
<path class="current" d="M805 285 L695 285"/>
<text class="body" x="713" y="307" fill="#d65a31">grid Ig,h</text>
<path class="current" d="M615 315 L585 340"/>
<text class="body" x="618" y="329" fill="#d65a31">shunt Ish,h</text>

<rect class="greenbox" x="70" y="515" width="860" height="72" rx="9"/>
<text class="head" x="93" y="544">Reported and reserved separately at every physical path</text>
<text class="formula" x="93" y="570">Isw,rms² = 2 Σh |Ih|²  →  converter, grid, shunt, and neutral RMS / peak-to-peak / total RMS</text>
</svg>"""
    write(path, svg)
end

function dc_source_sharing_figure(path)
    svg = """<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="600" viewBox="0 0 1000 600">
<rect width="1000" height="600" fill="#ffffff"/>
<style>
 .title { font:700 24px sans-serif; fill:#172033 }
 .sub { font:14px sans-serif; fill:#526174 }
 .head { font:700 16px sans-serif; fill:#172033 }
 .body { font:13px sans-serif; fill:#334155 }
 .formula { font:15px serif; fill:#172033 }
 .wire { stroke:#172033; stroke-width:2.5; fill:none }
 .blue { stroke:#1769aa; stroke-width:2.5; fill:none; marker-end:url(#blue) }
 .orange { stroke:#d65a31; stroke-width:2.5; fill:none; marker-end:url(#orange) }
 .green { stroke:#33865a; stroke-width:2.5; fill:none; marker-end:url(#green) }
 .box { fill:#f7f9fc; stroke:#8ca0b8; stroke-width:1.7 }
 .bluebox { fill:#edf7ff; stroke:#1769aa; stroke-width:1.8 }
 .warn { fill:#fff6eb; stroke:#d65a31; stroke-width:1.8 }
</style>
<defs>
 <marker id="blue" markerWidth="9" markerHeight="9" refX="8" refY="3" orient="auto"><path d="M0,0 L0,6 L9,3 z" fill="#1769aa"/></marker>
 <marker id="orange" markerWidth="9" markerHeight="9" refX="8" refY="3" orient="auto"><path d="M0,0 L0,6 L9,3 z" fill="#d65a31"/></marker>
 <marker id="green" markerWidth="9" markerHeight="9" refX="8" refY="3" orient="auto"><path d="M0,0 L0,6 L9,3 z" fill="#33865a"/></marker>
</defs>
<text class="title" x="42" y="40">DC carrier harmonics: source and capacitor current sharing</text>
<text class="sub" x="42" y="65">The bridge injects one correlated harmonic current; the physical parallel network decides where it closes</text>

<rect class="bluebox" x="55" y="165" width="190" height="135" rx="10"/>
<text class="head" x="80" y="197">shared-carrier bridge</text>
<text class="formula" x="82" y="230">Îbridge,h</text>
<text class="body" x="80" y="265">3-leg · split DC · 4-leg</text>
<path class="blue" d="M245 232 L365 232"/>
<circle cx="380" cy="232" r="7" fill="#172033"/>
<text class="body" x="342" y="212">DC node V̂h</text>

<path class="wire" d="M380 232 L380 120 L520 120"/>
<rect class="box" x="520" y="85" width="245" height="105" rx="9"/>
<text class="head" x="544" y="116">upstream source branch</text>
<text class="formula" x="544" y="148">Zs,h = Rs + jΩhLs</text>
<text class="body" x="544" y="174">optional; open by default</text>
<path class="wire" d="M765 120 L850 120 L850 395"/>
<path class="orange" d="M490 102 L425 102"/>
<text class="body" x="430" y="91" fill="#d65a31">Îsource,h</text>

<path class="wire" d="M380 232 L380 395 L520 395"/>
<rect class="box" x="520" y="330" width="245" height="130" rx="9"/>
<text class="head" x="544" y="361">DC-link capacitance</text>
<text class="formula" x="544" y="394">YC,h = jΩhCeq</text>
<text class="body" x="544" y="423">split: Ceq = CuCl/(Cu+Cl)</text>
<text class="body" x="544" y="446">V̂u=(Ceq/Cu)V̂ · V̂l=(Ceq/Cl)V̂</text>
<path class="wire" d="M765 395 L850 395"/>
<path class="green" d="M490 377 L425 377"/>
<text class="body" x="437" y="366" fill="#33865a">Îcap,h</text>

<rect class="warn" x="85" y="495" width="830" height="75" rx="9"/>
<text class="head" x="110" y="525">Harmonic KCL and antiresonance diagnostic</text>
<text class="formula" x="110" y="553">Îbridge,h + Îcap,h + Îsource,h = 0     ·     ρdc = min |YC,h+Ys,h|/(|YC,h|+|Ys,h|)</text>
</svg>"""
    write(path, svg)
end

mkpath(OUT)
waveform_figure(joinpath(OUT, "generated-ripple-waveforms.svg"))
utilisation_figure(joinpath(OUT, "generated-topology-utilisation.svg"))
lcl_circuit_figure(joinpath(OUT, "generated-lcl-circuit.svg"))
split_link_figure(joinpath(OUT, "generated-split-link-asymmetry.svg"))
pwm_ripple_figure(joinpath(OUT, "generated-pwm-ripple-closure.svg"))
ac_ripple_network_figure(joinpath(OUT, "generated-ac-ripple-network.svg"))
dc_source_sharing_figure(joinpath(OUT, "generated-dc-source-sharing.svg"))
println("Generated IBR figures in $OUT")
