predicate fzn_table_bool(array[int] of var bool: x,
                           array[int] of bool: t);

predicate fzn_table_bool(array[int] of var bool: x, array[int, int] of bool: t) =
  fzn_table_bool(x,array1d(t));

