## ONCVPSP CODING STANDARDS

The coding standards I've tried to follow and will want for any contributed additions are pretty basic.
They should also be obvious from examininga few of the source routines.

1) Fortran 95 free-format.

2) One subroutine per file.

3) All data is passed as subroutine scalar or array variables.  There are to
be no common blocks, data in modules, or derived data types.

4) All variables are to be explicitly typed (eg. "implicit none").

5) No one-letter variable names.

6) Input, output, and local variables should be grouped.

7) All local arrays should be allocated and deallocated (unless they're tiny).

8) Each subroutine should have a brief statement of what it does immediately
before or after the subroutine statement.

9) All subroutine arguments should be defined in comments following the
subroutine statement.

10) There should be comments throughout to indicate what is going on.  
Any new math should be described in a text or pdf contribution to the doc directory.

11) I'd like to stick to one input file read on standard input and one
output file to standard output, with shell scripts to extract any new
sections.

12) If you add new input data, the process of reading and checking it should
be done in such a way that older files still run.  This could also be
accomplished by building a separate executable for the new feature.

There are probably other conventions I've been following but you get the idea.
