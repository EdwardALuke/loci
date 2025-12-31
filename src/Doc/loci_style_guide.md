C++ Guidelines

* Put a space between the end of a line and the semicolon.

* Put spaces in comma-separated lists, like those in function arguments, function calls, and variable declarations.

* Use 80 characters as the maximum line length. If a line exceeds this length, break it before an operator.

* Indentation is 2 spaces.

* Structure blocks as follows:

void function() {
  if (condition) {
    // code
  } else {
    // code
  }
}

* If a c++ function has a very long argument list, indent by tabbing over (tab with 2 spaces) 3 times like shown below.
void function(argumentA, argumentB, argumentC, argumentD,
      argumentE, argumentF, argumentG) {


}

Loci Guidelines

* For long rules, break the rule before the 80 character limit, but do not indent the continuation line.

* Rule structures as follows (showing the line breaking style):

$rule pointwise(variableA<-variableB), constraint(constraintA),
options(disable_threading), {
  // code
}

* Use a /// comment on rules as they are parsed by the Loci-Pre-Processor. Extraneous comments in the source code should use //.
