
% 99 Prolog Problems
% Tyler Green's Solutions
% December 2007

%utilities

cons(X,Y,[X|Y]).

my_append([], L, L).
my_append([X|L1], L2, [X|L3]) :- my_append(L1, L2, L3).

foldr(_, Z, [], Z).	
foldr(P, Z, [X|Xs], R) :-
	foldr(P, Z, Xs, R1),
	Q =.. [P, X, R1, R],
	call(Q).

foldl(_, Z, [], Z).
foldl(P, Z, [X|Xs], R) :-
	Q =.. [P,Z,X,R1],
	call(Q),
	foldl(P, R1, Xs, R).

zip(_, [], []).
zip([], _, []).
zip([X|Xs], [Y|Ys], [[X|Y] | Z]) :- zip(Xs, Ys, Z).

take(0,_,[]).
take(N,[X|Y],[X|Xs]) :- N1 is N - 1, take(N1,Y,Xs).

% sorted tree dictionary

% ascii codes: newline = 10, space = 32
words(C,Ws) :- separate(32,C,Ws).
lines(C,Cs) :- separate(10,C,Cs).

separate(_,[],[]). 
separate(Sym,X,[L|R1]) :- break(Sym,X,L,R), separate(Sym,R,R1).
separate(Sym,X,[L]) :- break(Sym,X,L,[]). 

break(Sym,[],[],_).
break(Sym,[Sym|Xs],[],Xs).
break(Sym,[X|Xs],[X|L],R) :- break(Sym,Xs,L,R).

prefix(X,Y) :- append(Y,_,X).

islist([]).
islist([_|T]) :- islist(T).

length1([],0).
length1([_|T], N) :- length(T, N1), N1 is N - 1.

eq_length([],[]).
eq_length([_|Xs],[_|Ys]) :- eq_length(Xs, Ys). 

member(X, [X|_]).
member(X, [_|Ys]) :- member(X, Ys).

% if E is a new element, return a new tree with E in it.  
% otherwise E is already present, so fail
is_new_tree(E, leaf, tree(E, leaf, leaf)).
is_new_tree(E, tree(E1, L, R), tree(E1, L, R1)) :-
	E @< E1,
	is_new_tree(E, R, R1).
is_new_tree(E, tree(E1, L, R), tree(E1, L1, R)) :-
	E @> E1,
	is_new_tree(E, L, L1).

tree_mem(E, tree(E, _, _)) :- !.
tree_mem(E, tree(E1,L,_)) :-
	E @< E1,
	tree_mem(E, L).
tree_mem(E, tree(E1,_,R)) :-
	E @> E1,
	tree_mem(E,R).


replace(X, Y, [], []).
replace(X, Y, [X|T], [Y|T]) :- !.
replace(X, Y, [E|T], [E|R]) :- replace(X, Y, T, R).

replace1(X, Y, [], _, []).
replace1(X, Y, [X|T], [_|T1], [Y|T1]).
replace(X, Y, [E|T], [E1|T1], [E1|R]) :- replace(X, Y, T, T1, R).

switch(X, Y, X, Y) :- !.
switch(X, Y, Z, Z).

alter(X, Y, [], [], []).
alter(X, Y, [X|Es], [_|Zs], [Y|Zs]).
alter(X, Y, [E|Es], [Z|Zs], [Z|R]) :- \+ X = E, alter(X, Y, Es, Zs, R).

map(_, [], []).
map(P, [X|L], [Y|M]) :-
	Q =.. [P,X,Y], call(Q), map(P, L, M).

map(P, [], _, []).
map(P, _, [], []).
map(P, [X|Xs], [Y|Ys], [Z|Zs]) :-
	Q =.. [P,X,Y,Z],
	call(Q),
	map(P, Xs, Ys, Zs).

% what if this backtracks?
map1(Pred, []).
map1(Pred, [X|Xs]) :- Q =.. [Pred, X], call(Q), map1(Pred, Xs).

applist(_, []).
applist(P, [X|Xs]) :- Q =.. [P, X], call(Q), applist(P,Xs).

fmap(X, Y, Z) :- map(X, Y, Z).
fmap(X, Y) :- applist(X, Y).


mult(X, Y, Z) :- Z is X * Y.

fac(X,Y) :- range(1, X, R), foldr(mult, 1, R, Y).

%filter  - use include, exclude
filter(X,Y,Z) :- include(X,Y,Z).

pprint(X) :- fmap(put, X), nl.

%%%% String utilties

wspaces(0) :- !.
wspaces(N) :- write(' '), N1 is N - 1, wspaces(N1). 

% 0 all but last element of the list

butLast([X], []).
butLast([X|Xs], [X|R]) :- butLast(Xs, R).

% 1 last element of a list

last1([X],X).
last1([_|T], X) :- last1(T, X).

% 2 second to last element of a lit

second_last(X, [X, _]).
second_last(X, [_|T]) :- second_last(X, T).

% 3 find the k'th element of a list

kth(0, [X|_], X) :- !.
kth(N, [_|T], X) :- N1 is N - 1, kth(N1, T, X).

% 4 count the number of elements in a list
count([], 0).
count([_|T], N) :- count(T, N1), N is N1 + 1. 

						% 5 reverse a list

reverse(L, X) :- revaux(L, [], X).
rev(L, X) :- revaux(L, [], X).

revaux([], A, A).
revaux([X|Xs], A, N) :- revaux(Xs, [X|A], N).

% 6 is X a palindrome?
is_palindrome(X) :- rev(X,X).


flattenT(leaf, []).
flattenT(tree(E, L, R), List) :- 
    flattenT(L, L1),
    flattenT(R, R1),
    append(L1, [E|R1], List).

% 7 flatten(input, output).
flatten([], []).
flatten([X|XS], Z) :- atom(X), flatten(XS, Z1), append([X], Z1, Z).
flatten([X|XS], Z) :- flatten(X, Z1), flatten(XS, Z2), append(Z1, Z2, Z).

% this doesn't seem to work
flatten1([], []).
flatten1([X|XS], Z) :- atom(X), append([X], flatten(XS, _), Z).
flatten1([X|XS], Z) :- append(flatten(X, _), flatten(XS, _), Z).

% 8 Eliminate consecutive duplicates of list elements.

compress([],[]).
compress([X,X|XS], Z) :- compress([X|XS], Z).
compress([X|XS], Z) :- compress(XS, Z1), append([X], Z1, Z).

% pack consecutive duplicates into lists

pack([],[]).
pack([X|Xs],[Z|Zs]) :- transfer(X,Xs,Ys,Z), pack(Ys,Zs).

transfer(X,[],[],[X]).
transfer(X,[Y|Ys],[Y|Ys],[X]) :- X \= Y.
transfer(X,[X|Xs],Ys,[X|Zs]) :- transfer(X, Xs, Ys, Zs).

% 10 run-length encoding

encode(L, X) :- pack(L, PL), map_count(PL, X).

map_count([],[]).
map_count([X|Xs], [Y|Ys]) :- count(X,C), [M|Ms]=X, Y=[M,C], map_count(Xs,Ys).

% 12 Encode/Decode a run-length encoding
% this can both encode and decode

decode([],[]).
decode([[Sym,N]|Xs], [Y|Ys]) :- n_dup(Sym,N,Y), decode(Xs,Ys).

code1(Encoded, Decoded) :- maplist(n_dup, Encoded, Decoded).

% duplicates X (N times)
n_dup(X, 0, []).
n_dup(X, N, [X|Xs]) :- n_dup(X, N1, Xs), N is N1 + 1.

% mapping
map([],[]).
map([X|Xs], [Y|Ys]) :- fn(X, Y), map(Xs, Ys).

% 14 duplicate each elem. of a list
dup([],[]).
dup([X|XS], [X,X|YS]) :- dup(XS, YS).

double(X, [X,X]).

% 16 Drop the first k elements

drop(_, [], _) :- error( 'too many to drop').
drop(0, X, X).
drop(N, [X|Xs], Z) :- N1 is N -1, drop(N1,Xs,Z).

%drop([X|Xs], 0, N, [Y|Ys]) :- 
%drop([X|Xs], I, N, [Y|Ys]) :- drop1(

drop1(_, 0, []).
%drop1([X|Xs], N, [Y|Ys]) :- drop(

% 17 split a list a position k
% A more efficient solution is definitely possible

split(L, K, L1, L2) :- append(L1, L2, L), length(L1, K).

split1(L,0,[],L).
split1([X|Xs],K,[X|Ys],Zs) :- K1 is K - 1, split1(Xs, K1, Ys, Zs).

% 18 extract a slice from a list

slice1(I, K, L, X) :-
	I1 is I -1,
	split(L,I1,_,L1),
	N is K - I,
	split(L1,N,X,_).


%% includes the S0th element, but not the SNth element => (S0, SN]
slice(S0, SN, L, L1) :-
	append(Prefix,Rhs,L),
	length(Prefix,S0),
	append(L1,Suffix,Rhs),
	S1 is SN - S0,
	length(L1, S1).

% 19 Rotate a list N places to the left

rotate(L, N, X) :- split(L, N, L1, L2), append(L2, L1, X).

% 20 remove the k'th element from the list

remove_at(X,[X|Xs],0, Xs).
remove_at(X,[Y|Xs],K, [Y|Ys]) :- K > 0, K1 is K - 1, remove_at(X,Xs, K1, Ys).

% 21 insert an element in to a list at position k

nthcdr([X|XS], 0, XS).
nthcdr([X|XS], N, Z) :- nthcdr(XS, N1, Z), N is N1 + 1.

num(0).
num(s(X)) :- num(X).

sum(X, 0, X) :- num(X).
sum(X, s(Y), s(Z)) :- sum(X, Y, Z).

% 22 create a list containing all the integers within a given range

range(N,N, [N]).
range(N1, N2, [N1|Ns]) :- N3 is N1 + 1, range(N3, N2, Ns).

% 23
%use_module(library(random)).

rnd_select(_,0,[]).
rnd_select(List, N, [X|Xs]) :- 
	count(List, Length),
	random(0, Length, K),
	remove_at(X,List,K,L2),
	N1 is N - 1, rnd_select(L2, N1, Xs).

% 26 N Choose K Combinations

% 31 Determine whether a given number is prime

% 33

% binary tree? 

istree(nil).
istree(t(_,L,R)) :- istree(L), istree(R).

% My tree representation
istree1([]).
istree1([A|B]) :- istree1(A), istree1(B).

% count the leaves of a tree

% repeatedly execute Goal N times (good for profiling small preds)
repeat(_, 0) :- !.
repeat(Goal, N) :-
	call(Goal),
	N1 is N - 1,
	repeat(Goal, N1).


% Chp 3.6 Programming in Prolog

assembly(bike, [wheel,wheel,frame]).
assembly(wheel, [spoke,rim,hub]).
assembly(frame, [rearframe, frontframe]).
assembly(frontframe, [fork, handles]).
assembly(hub, [gears, axle]).
assembly(axle, [bolt, nut]).

basicpart(rim).
basicpart(spoke).
basicpart(rearframe).
basicpart(handles).
basicpart(gears).
basicpart(bolt).
basicpart(nut).
basicpart(fork).

partsof(X,[X]) :- basicpart(X).

partsof(X,P) :-
	assembly(X, Subparts),
	partsoflist(Subparts, P).

partsoflist([],[]).
partsoflist([P|Tail], Total) :-
	partsof(P, Headparts),
	partsof(Tail, Tailparts),
	append(Headparts, Tailparts, Total).

% scheme, haskell, prolog, ruby, forth, 

% lisp family (chicken, plt, arc)
% functional (haskell ocaml)
% stack (factor, forth, joy)
% oo (smalltalk)
% imperative (ruby, c, python)





%% Local Variables: 
%% mode: CIAO
%% update-version-comments: "off"
%% End:

