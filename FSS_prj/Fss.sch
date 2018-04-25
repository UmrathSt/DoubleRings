<Qucs Schematic 0.0.18>
<Properties>
  <View=90,40,1387,1102,0.937879,0,0>
  <Grid=10,10,1>
  <DataSet=Fss.dat>
  <DataDisplay=Fss.dpl>
  <OpenDisplay=1>
  <Script=Fss.m>
  <RunScript=0>
  <showFrame=0>
  <FrameText0=Title>
  <FrameText1=Drawn By:>
  <FrameText2=Date:>
  <FrameText3=Revision:>
</Properties>
<Symbol>
</Symbol>
<Components>
  <C C1 1 490 260 17 -26 0 1 "1 pF" 1 "" 0 "neutral" 0>
  <L L2 1 490 360 10 -26 0 1 "0.498 nH" 1 "" 0>
  <GND * 1 490 390 0 0 0 0>
  <GND * 1 360 390 0 0 0 0>
  <GND * 1 740 390 0 0 0 0>
  <Eqn absS 1 960 120 -26 16 0 0 "Smeta=abs(S[1,1])" 1 "yes" 1>
  <TLIN Line1 1 740 260 20 -26 0 1 "179.7 Ohm" 1 "1 m" 1 "0 dB" 0 "26.85" 0>
  <Eqn Z 1 970 250 -26 16 0 0 "Absorption=1-abs(S[1,1])^2" 1 "yes" 0>
  <R R2 1 490 170 15 -26 0 1 "0.673 Ohm" 1 "26.85" 0 "0.0" 0 "0.0" 0 "26.85" 0 "european" 0>
  <Pac P1 1 360 290 18 -26 0 1 "1" 1 "376.7 Ohm" 1 "0 dBm" 0 "5 GHz" 0 "26.85" 0>
  <.SP SP1 1 300 650 0 78 0 0 "lin" 1 "2 GHz" 1 "8 GHz" 1 "10000" 1 "no" 0 "1" 0 "2" 0 "no" 0 "no" 0>
</Components>
<Wires>
  <360 120 490 120 "" 0 0 0 "">
  <490 290 490 330 "" 0 0 0 "">
  <360 120 360 260 "" 0 0 0 "">
  <360 320 360 390 "" 0 0 0 "">
  <490 120 490 140 "" 0 0 0 "">
  <490 200 490 230 "" 0 0 0 "">
  <490 120 740 120 "" 0 0 0 "">
  <740 120 740 230 "" 0 0 0 "">
  <740 290 740 390 "" 0 0 0 "">
</Wires>
<Diagrams>
  <Rect 832 790 488 407 3 #c0c0c0 1 00 1 2e+09 2e+09 8e+09 1 -0.0998927 0.5 1.09994 1 -1 1 1 315 0 225 "" "" "">
	<"Absorption" #ff0000 0 3 0 0 0>
  </Rect>
</Diagrams>
<Paintings>
</Paintings>
