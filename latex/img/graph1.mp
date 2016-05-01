prologues:=3;
verbatimtex
%&latex
\documentclass{minimal}
\begin{document}
etex
beginfig(0);

u = 2.3cm;
pair U,V,W,Z, sft;
path circle; 
circle = fullcircle scaled .2u;

U = (0,0);
V = (u,0);
W = (u,u);
Z = (0,u);
sft = (0,.1u);

fill circle withcolor .8*white;



draw U--Z--W--V--Z;
draw U withpen pencircle scaled .19u withcolor white;
draw W withpen pencircle scaled .19u withcolor white;
draw V withpen pencircle scaled .19u withcolor white;
draw Z withpen pencircle scaled .19u withcolor white;
draw circle shifted U;
draw circle shifted V;
draw circle shifted W;
draw circle shifted Z;

label.bot(btex 1 etex,U+sft);
label.bot(btex 2 etex,V+sft);
label.bot(btex 3 etex,W+sft);
label.bot(btex 4 etex,Z+sft);

label.lft(btex 4 etex,1/2[U,Z]);
label.rt(btex 4 etex,1/2[V,W]);
label.top(btex 4 etex,1/2[W,Z]);
label.urt(btex 8 etex,1/2[V,Z]);
endfig;
end