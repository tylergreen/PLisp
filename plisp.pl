%% PLisp, a scheme interpreter in prolog
%% Tyler Green
%% Summer 2008 continued

%% This is the basic scheme interpreter.
%% Lacks support in the following: strings, continuations
%% Use this as the underlying foundation of more advanced versions.

%% run_scheme/2 evaluates a single scheme expression represented as prolog lists.
%% example run_scheme([cons, 2, [quote, [4,6]]], R).   R = [2, 4, 6]

%% extensions:  improve pretty lisp printer. spacing is ugly
%%              better error checking/ failure other than just "no"
%%              implement primitive procedure better so can use '+' instead of 'plus'
%%              tail recursion
%%              first class continuations
%%              efficient implementation
%%              R5R6 or R6R6 compliance

%% Notes remember your syntax you implemented, watch for bugs in the primitve procedures, 
% especially the aritmetical ones

eval(Exp, Env, Exp, Env) :- self_eval(Exp).                            %done
eval(Exp, Env, V, Env) :- variable(Exp), lookup_var_val(Exp, Env, V). %done
eval(Exp, Env, V, Env) :- quoted(Exp,V).                        %done      
eval(Exp, Env, V, Env1) :- eval_assignment(Exp, Env, V, Env1).  %done    
eval(Exp, Env, V, Env1) :- eval_defn(Exp, Env, V, Env1).        %done    
eval(Exp, Env, V, Env1) :- eval_if(Exp, Env, V, Env1).          %done
eval(Exp, Env, V, Env)  :- make_proc(Exp, Env, V).              %done        
eval(Exp, Env, V, Env1) :- eval_sequence(Exp, Env, V, Env1).    %done       
eval(Exp, Env, V, Env1) :- eval_cond(Exp, Env, V, Env1).        %done             
eval(Exp, Env, V, Env1) :- eval_app(Exp, Env, V, Env1).         %done             

%%% top level

plisp :- init_env(Env), repl(Env).

repl(Env) :-
	write('p-lisp> '), 
	lisp_read(X),
	eval(X, Env, V, Env1),
	lisp_print(V),
	repl(Env1).

lisp_read(X) :- read(X).
lisp_print(X) :- writeln(X).

%%%%% self evaluating atoms
self_eval(nil).
self_eval(Exp) :- number(Exp).

%%% variables

variable(X) :-  atom(X).

%%%% quotation

quoted([quote, T], T).

%%%% definition

eval_defn([define, Var, Expr], Env, Var, Env2) :-
	!,
	eval(Expr, Env, Value, Env1),
	define_var(Var, Value, Env1, Env2).

%%% assignment isn't working properly !! 

eval_assignment([set, Var, Val], Env, Var, Env2) :-
	!,
	eval(Val, Env, Value, Env1),
	set_var_val(Var, Value, Env1, Env2).

%%%% if special form
eval_if([if, Pred, Conseq, Alter], Env, R, E2) :-
	eval(Pred, Env, R1, E1),
	!, truep(R1),
	eval(Conseq, E1, R, E2).
eval_if([if, Pred, Conseq, Alter], Env, R, E1) :-
	!, eval(Alter, Env, R, E1).

falsep(nil).
truep(Exp) :- \+ falsep(Exp).

%%%% apply
%% Important lisp note:  When the args get evaluated in a closure,
%% they are evaluated inside the closures environment.
eval_app([Operator | Operands], Env, R, Env1) :- 
	eval(Operator, Env, Proc, _), 
	list_of_vals(Operands, Env, Args, Env1),
	apply(Proc, Args, R).

%structure of a proc:  procedure(Params, Body, Env) 

apply(compound_proc(Params, Body, Env), Args, R) :-
	!,
	extend_env(Params, Args, Env, Env1),
	eval_seq(Body, Env1, R, _). 
apply(Proc, Args, R) :- apply_primitive_proc(Proc, Args, R).

%%primitive_proc(Proc) :-
%% List of prim procs = [car, cdr, cons, list, append, null, plus, minus, product, foldr].

apply_primitive_proc(Proc, Args, R) :- P =.. [Proc, Args, R], call(P).

%%% cond 
%(cond ((= a b) c)
%      (else (- a d)))

%[cond, [[=,a,b], c],
%	[else, [-, a, d]]]


eval_cond([cond|Body], Env, R, Env1) :-
	cond_if(Body, Nested_Ifs),
	eval(Nested_Ifs, Env, R, Env1).

% convert cond statement to nested if statements
cond_if([], nil).  %no else clause
cond_if([[else|Actions]], R) :- !, sequence_exp(Actions, R). %error checking opportunity here
cond_if([C|Cs], [if, Pred, Conseq, Altern]) :- 
	[Pred|Actions] = C,
	sequence_exp(Actions, Conseq),
	cond_if(Cs, Altern).

% add implicit begin to 'then' part of cond clause
sequence_exp([],[]).
sequence_exp([Exp], Exp).
sequence_exp(Exprs, [begin|Exprs]).

% the recursion can be abstracted into map/4. 
list_of_vals([], Env, [], Env).
list_of_vals([E|Exps], Env, [R|Rs], Env2) :-
	eval(E, Env, R, Env1),
	list_of_vals(Exps, Env1, Rs, Env2).

% =.. approach?
%list_of_vals(Exp, Env, R) :- Eval =.. [eval, Env, maplist(Eval, Env, 

%%%% lambda
% (\ (x) (list (+ 5 x)))

compound_proc(Params, Body, Env).
make_proc([\, Params|Body], Env, compound_proc(Params, Body, Env)).

%%%% primitive procedures (scheme procedures written in prolog) 
%% remember takes a single list as arg

car([H|_], H).
cdr([_|T], T).
cons([X,Xs], [X|Xs]).
list(L,L).

% want to do arbitrary length append.
append([[], L], L).
append([[X|L1], L2], [X|L3]) :- append(L1, L2, L3).

null([]).

%% add more
eq(X, X).
listp([]).
listp([X|_]).

plus(Ns, Sum) :- foldr(bin_plus, 0, Ns, Sum).
minus([S|Ss], Diff) :- plus(Ss, Sum), Diff is S - Sum. 
product(Ns, Product) :- foldr(bin_mult, 1, Ns, Product).

foldr(Op, Init, [], Init).
foldr(Op, Init, [X|Xs], R) :- foldr(Op, Init, Xs, R1), P =.. [Op, X, R1, R], call(P).

bin_plus(A,B, Sum) :- Sum is A + B.
bin_mult(A,B, Prod) :- Prod is A * B.

%%%% Sequences
eval_sequence([begin|Exprs], Env, R, Env1) :- eval_seq(Exprs, Env, R, Env1).

eval_seq([Exp], Env, R, Env1) :- eval(Exp, Env, R, Env1).
eval_seq([E|Exps], Env, R, NewEnv) :-  
	eval(E, Env, _, Env1),          %throw out value
	eval_sequence(Exps, Env1, R, NewEnv).

%%% Predicates on Environments
% an environment is a list of frames  Env :: [Frame]
% a frame is a pair of lists:  
% returns the val bound to symbol Var in the environment Env, signals an error? if unbound
 
env([F|Frames]).
frame(Vars, Vals, [Vars|Vals]).

enclosing_env([E|Es], Es).
is_empty_env([]).

init_env(X) :- 
	L = [car, cdr, cons, list, append, null, plus, minus, product, foldr],
	is_empty_env(Empty),
	extend_env(L, L, Empty, X).

lookup_var_val(Var, Env, Val) :- 
	first_frame(Env, F),
	frame(Vars, Vals, F),
	scan(Var, Vars, Vals, Env, Val).

set_var_val(Var, Val, Old_Env, New_Env) :-
	update_env(Var, Val, Old_Env, New_Env).

define_var(Var, Val, [F|Fs], [F1|Fs]) :-
	add_binding_to_frame(Var, Val, F, F1).

extend_env(Vars, Vals, Base_Env, [F|Base_Env]) :-
	eq_length(Vars, Vals),
	frame(Vars, Vals, F).
extend_env(_,_,_,_) :- throw('incorrect number of args supplied').

first_frame([F|_], F).

add_binding_to_frame(Var, Val, [Vars|Values], [[Var|Vars]|[Val|Values]]).

scan(Var, [], [], Env, R) :- 
	enclosing_env(Env, Enc_Env),
	lookup_var_val(Var, Enc_Env, R).
scan(Var, [Var|_], [Val|_], _, Val).
scan(Var, [_|Vars], [_|Vals], Env, R) :- scan(Var, Vars, Vals, Env, R).

update_env(Var, Val, [F|Fs], [F1|Fs]) :-
	frame(Vars, Vals, F),   %frame won't backtrack, fails after first
	member(Var, Vars),
	update_frame(Var, Val, Vars, Vals, New_Vals),
	frame(Vars, New_Vals, F1).
% if Var is not present in this frame go, to the next frame (the enclosing env)
update_env(Var, Val, [F|Fs], [F|New_Env]) :-  
	update_env(Var, Val, Fs, New_Env).
update_env(_,_,[],_) :-
	throw('unbound variable').

% update frame is true if X is present in the Vars and Y is the new value for X in Es
update_frame(Var, NewVal, [], [], []).
update_frame(Var, NewVal, [Var|Vars], [_|Vals], [NewVal|Vals]).
update_frame(Var, NewVal, [X|Vars], [Y|Vals], [Y|R]) :-
	Var \= X,
	update_frame(Var, Val, Vars, Vals, R).


