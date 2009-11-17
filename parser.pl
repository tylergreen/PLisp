% Lisp grammar using Definite Clause Grammars

read_sexp(S) --> spaces, sexp(S), spaces.

sexp(S) --> list(S). 
sexp(S) --> atom(S).

list(Ns) --> lparen, spaces, nodes(Ns), spaces, rparen.

nodes([N|Ns]) --> sexp(N), spaces, nodes(Ns).
nodes([]) --> [].

lparen --> [40].
rparen --> [41].

spaces --> [32], spaces.
spaces --> [].

atom(A) --> integer(A).

integer(I) -->
        digit(D0),
        digits(D),
        { number_chars(I, [D0|D])
        }.

digits([D|T]) -->
        digit(D), !,
        digits(T).
digits([]) -->
        [].

digit(D) -->
        [D],
        { code_type(D, digit)
        }.




	
	
	
