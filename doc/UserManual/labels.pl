# LaTeX2HTML 2002-2-1 (1.71)
# Associate labels original text with physical files.


$key = q/intro/;
$external_labels{$key} = "$URL/" . q|node2.html|; 
$noresave{$key} = "$nosave";

$key = q/tutorial/;
$external_labels{$key} = "$URL/" . q|node3.html|; 
$noresave{$key} = "$nosave";

$key = q/cite_Giacovazzo92/;
$external_labels{$key} = "$URL/" . q|node7.html|; 
$noresave{$key} = "$nosave";

$key = q/ewald-construction/;
$external_labels{$key} = "$URL/" . q|node6.html|; 
$noresave{$key} = "$nosave";

1;


# LaTeX2HTML 2002-2-1 (1.71)
# labels from external_latex_labels array.


$key = q/intro/;
$external_latex_labels{$key} = q|1|; 
$noresave{$key} = "$nosave";

$key = q/tutorial/;
$external_latex_labels{$key} = q|1.1|; 
$noresave{$key} = "$nosave";

$key = q/ewald-construction/;
$external_latex_labels{$key} = q|3|; 
$noresave{$key} = "$nosave";

1;

