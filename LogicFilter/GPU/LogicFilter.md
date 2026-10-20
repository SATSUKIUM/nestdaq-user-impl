# LogicFilter

GPUの処理の流れ・データ構造・同期については[実装説明](Implementation.md)を参照。

## Overview
LogicFilter processes input TDC (Time-to-Digital Converter) data based on a user-defined logical expression. If the expression evaluates to "true" at a specific timestamp, that timestamp is added to the output data.
The processing time unit is four times the unit of the input TDC data. For example, if the input TDC data unit is 1 ns, the output data unit will be 4 ns.

## Redis DB Keys

### Key prefix
- Defaut: parameter:LogicFilter:

### trigger-signals
Defines inputs as `(group_id module_IP channel offset)` or
`(group_id-subgroup_id module_IP channel offset)`.
Group IDs are 0–31. Each subgroup contains up to 32 signals; inputs in a subgroup
are combined with AND, and subgroups assigned to the same group are combined with OR.
This allows more than 32 input signals in total.
- Default: `(0 0xc0a802a9 0 0) (1 0xc0a802a9 1 0)`

### trigger-expression
Defines the logical expression applied to the input signals. It supports both Infix notation and RPN (Reverse Polish Notation).

**Operands:**
- Numbers represent group IDs.

**Operators:**
- & : Logical AND
- | : Logical OR
- ! : Logical NOT

**Syntax:**
- Infix: Use parentheses ( and ) for order of operations.
- RPN: Must be prefixed with the keyword RPN.

- Default: RPN 0 1 &  
(Generates a trigger based on the coincidence of signal id 0 and 1.)

### trigger-width
Sets the TDC coincidence window. Let $T$ be the trigger-width value and $t$ be the TDC value; the window is defined as: 
$$[t - T/2, t + T/2]$$
The time unit for $T$ is four times the time unit of the input TDC data.
- Default: 10
