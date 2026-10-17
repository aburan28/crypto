---- MODULE Recombination_TTrace_1792115830 ----
EXTENDS Sequences, Recombination, TLCExt, Toolbox, Naturals, TLC

_expression ==
    LET Recombination_TEExpression == INSTANCE Recombination_TEExpression
    IN Recombination_TEExpression!expression
----

_trace ==
    LET Recombination_TETrace == INSTANCE Recombination_TETrace
    IN Recombination_TETrace!trace
----

_inv ==
    ~(
        TLCGet("level") = Len(_TETrace)
        /\
        phase = ("Giant")
        /\
        allowed = ({0})
        /\
        count = (2)
        /\
        j = (2)
        /\
        stride = (1)
        /\
        point = (0)
        /\
        n = (1)
        /\
        hits = ({0})
        /\
        result = (-1)
        /\
        truth = (0)
        /\
        round = (0)
        /\
        survivors = (0..1)
        /\
        baby = ((0 :> {0}))
        /\
        giant = (1)
    )
----

_init ==
    /\ phase = _TETrace[1].phase
    /\ stride = _TETrace[1].stride
    /\ point = _TETrace[1].point
    /\ truth = _TETrace[1].truth
    /\ j = _TETrace[1].j
    /\ n = _TETrace[1].n
    /\ baby = _TETrace[1].baby
    /\ survivors = _TETrace[1].survivors
    /\ hits = _TETrace[1].hits
    /\ count = _TETrace[1].count
    /\ allowed = _TETrace[1].allowed
    /\ result = _TETrace[1].result
    /\ round = _TETrace[1].round
    /\ giant = _TETrace[1].giant
----

_next ==
    /\ \E i,j \in DOMAIN _TETrace:
        /\ \/ /\ j = i + 1
              /\ i = TLCGet("level")
        /\ phase  = _TETrace[i].phase
        /\ phase' = _TETrace[j].phase
        /\ stride  = _TETrace[i].stride
        /\ stride' = _TETrace[j].stride
        /\ point  = _TETrace[i].point
        /\ point' = _TETrace[j].point
        /\ truth  = _TETrace[i].truth
        /\ truth' = _TETrace[j].truth
        /\ j  = _TETrace[i].j
        /\ j' = _TETrace[j].j
        /\ n  = _TETrace[i].n
        /\ n' = _TETrace[j].n
        /\ baby  = _TETrace[i].baby
        /\ baby' = _TETrace[j].baby
        /\ survivors  = _TETrace[i].survivors
        /\ survivors' = _TETrace[j].survivors
        /\ hits  = _TETrace[i].hits
        /\ hits' = _TETrace[j].hits
        /\ count  = _TETrace[i].count
        /\ count' = _TETrace[j].count
        /\ allowed  = _TETrace[i].allowed
        /\ allowed' = _TETrace[j].allowed
        /\ result  = _TETrace[i].result
        /\ result' = _TETrace[j].result
        /\ round  = _TETrace[i].round
        /\ round' = _TETrace[j].round
        /\ giant  = _TETrace[i].giant
        /\ giant' = _TETrace[j].giant

\* Uncomment the ASSUME below to write the states of the error trace
\* to the given file in Json format. Note that you can pass any tuple
\* to `JsonSerialize`. For example, a sub-sequence of _TETrace.
    \* ASSUME
    \*     LET J == INSTANCE Json
    \*         IN J!JsonSerialize("Recombination_TTrace_1792115830.json", _TETrace)

=============================================================================

 Note that you can extract this module `Recombination_TEExpression`
  to a dedicated file to reuse `expression` (the module in the 
  dedicated `Recombination_TEExpression.tla` file takes precedence 
  over the module `Recombination_TEExpression` below).

---- MODULE Recombination_TEExpression ----
EXTENDS Sequences, Recombination, TLCExt, Toolbox, Naturals, TLC

expression == 
    [
        \* To hide variables of the `Recombination` spec from the error trace,
        \* remove the variables below.  The trace will be written in the order
        \* of the fields of this record.
        phase |-> phase
        ,stride |-> stride
        ,point |-> point
        ,truth |-> truth
        ,j |-> j
        ,n |-> n
        ,baby |-> baby
        ,survivors |-> survivors
        ,hits |-> hits
        ,count |-> count
        ,allowed |-> allowed
        ,result |-> result
        ,round |-> round
        ,giant |-> giant
        
        \* Put additional constant-, state-, and action-level expressions here:
        \* ,_stateNumber |-> _TEPosition
        \* ,_phaseUnchanged |-> phase = phase'
        
        \* Format the `phase` variable as Json value.
        \* ,_phaseJson |->
        \*     LET J == INSTANCE Json
        \*     IN J!ToJson(phase)
        
        \* Lastly, you may build expressions over arbitrary sets of states by
        \* leveraging the _TETrace operator.  For example, this is how to
        \* count the number of times a spec variable changed up to the current
        \* state in the trace.
        \* ,_phaseModCount |->
        \*     LET F[s \in DOMAIN _TETrace] ==
        \*         IF s = 1 THEN 0
        \*         ELSE IF _TETrace[s].phase # _TETrace[s-1].phase
        \*             THEN 1 + F[s-1] ELSE F[s-1]
        \*     IN F[_TEPosition - 1]
    ]

=============================================================================



Parsing and semantic processing can take forever if the trace below is long.
 In this case, it is advised to uncomment the module below to deserialize the
 trace from a generated binary file.

\*
\*---- MODULE Recombination_TETrace ----
\*EXTENDS IOUtils, Recombination, TLC
\*
\*trace == IODeserialize("Recombination_TTrace_1792115830.bin", TRUE)
\*
\*=============================================================================
\*

---- MODULE Recombination_TETrace ----
EXTENDS Recombination, TLC

trace == 
    <<
    ([phase |-> "Baby",allowed |-> {0},count |-> 2,j |-> 0,stride |-> 1,point |-> 0,n |-> 1,hits |-> {},result |-> -1,truth |-> 0,round |-> 0,survivors |-> 0..1,baby |-> (0 :> {}),giant |-> 0]),
    ([phase |-> "Baby",allowed |-> {0},count |-> 2,j |-> 1,stride |-> 1,point |-> 0,n |-> 1,hits |-> {},result |-> -1,truth |-> 0,round |-> 0,survivors |-> 0..1,baby |-> (0 :> {0}),giant |-> 0]),
    ([phase |-> "Baby",allowed |-> {0},count |-> 2,j |-> 2,stride |-> 1,point |-> 0,n |-> 1,hits |-> {},result |-> -1,truth |-> 0,round |-> 0,survivors |-> 0..1,baby |-> (0 :> {0}),giant |-> 0]),
    ([phase |-> "Giant",allowed |-> {0},count |-> 2,j |-> 2,stride |-> 1,point |-> 0,n |-> 1,hits |-> {},result |-> -1,truth |-> 0,round |-> 0,survivors |-> 0..1,baby |-> (0 :> {0}),giant |-> 0]),
    ([phase |-> "Giant",allowed |-> {0},count |-> 2,j |-> 2,stride |-> 1,point |-> 0,n |-> 1,hits |-> {0},result |-> -1,truth |-> 0,round |-> 0,survivors |-> 0..1,baby |-> (0 :> {0}),giant |-> 1])
    >>
----


=============================================================================

---- CONFIG Recombination_TTrace_1792115830 ----
CONSTANTS
    MaxGroup = 3
    MaxCount = 5
    MaxRounds = 1
    Legacy = TRUE

INVARIANT
    _inv

CHECK_DEADLOCK
    \* CHECK_DEADLOCK off because of PROPERTY or INVARIANT above.
    FALSE

INIT
    _init

NEXT
    _next

CONSTANT
    _TETrace <- _trace

ALIAS
    _expression
=============================================================================
\* Generated on Fri Oct 16 01:57:10 UTC 2026