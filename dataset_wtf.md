Just glancing at the documentation, you get the strong impression that
the 'dataset` class in the Statistics toolbox is a shoddy and
poorly-thought-out ripoff of R's `data.frame` class. This is true, and
true to form, Mathworks seriously fucked up several things about it.

Example:

###Grouping operation

One of the is a grouping operation. You have a dataset that veries
among two or three categirical axes (or binned axes) and you want to
plit it up, do stuff to each piece, and combine the results. you know,
the split-apply-combine operation (cf. WIckham) that is ubiquitous in
all data analysis. The only thing in the toolbox that approaches that
use case for datasets is `grpstats`, and, glancing at the help for
`grpstats` you see that you can use arbitrary functions, so you think
for a brief moment that using the `dataset` of MATLAB won't be so bad.

Not so fast. The only kind of operation that --- so grpstats splits of
by group, like you want, then splits up by column, which no one
wants. When have you ever wanted to the same statistics for all
columns in a dataset?! Different columns, generally, mean different
things! Some columns are experimenter-controlled variables, some
columns are observed data, some columns are categorical, some are
continuous. The situation of wanting to apply the same statistic to
all columns of a datatset arises approximately _never_ outside of
textbook examples.

So, actually, `grpstats` gets you nowhere along the way to doing what
you want with your data, and if you want womething like R's `tapply`
or even better `ddply` from the wonderful `plyr` package you have to
code it up yourself.

###Joins

The default join produced by `dataset/join(A,B)` is not an inner join,
or rather, it's an inner join the arbitrarily rejects cases; volating
the principle of least surprise, it's an _asymmetric_ operation; it
insists that all keys in B are matched my some key in A, while also
insisting that all keys of B have unique values.

Let's just look through sqlzoo.com, and see how often the default join would

So, let's set up one of hte simplest join operations imaginable, you
have a table of names `firstname, lastname, personID, householdID` and
an array of `householdID,

What are some of the operations you would like to do with this? Well,
let's say you have someone's name and you want to look their address.

[....expand on this example]

And if the intent is to make sure that all values are matched, it
doesn't work. Consider: (3 out of 4 in a matrix)

[....expand this example]

So, who the fuck knows what the default behavior of matlab's `join` is
good for. To get the join that everyone who's glanced at a database in
their life expects, you have to use `join(..., 'Type', 'inner')`.

But that's not all. When you "join" two datasets, if two non-key
values between the datasets overlap, it will rename them `key_left`
and `key_right`.  Fair enough; it has to do something. But in an inner
join, _by definition_ both left and right key values match. And the
default behavior is to use all matching field names as keys. But
MATLAB stull duplicates the column and renames it It `key_left` and
`key_right` for no reason even though both columns contain the same fucking
values _by definition_.

So, when you try using `@dataset/join` to do what normal people think
of as a 'join', 90% of the time you want to tack on `'Type', 'inner',
'MergeKeys', true` to your join arguments,, in the same fucking stupid
way that 90% of the time you `cellfun` you need to tack on a
`'UniformOutput', 0`.

I mean, if Mathworks wanted to rip off `data.frame` or other things from
R, you'd think that they would at least _try_ to translate some R code
into MATLAB just to see how many times more code it takes to do it in
MATLAB? They might have come up with some better default behaviors.


Also, check out this insane dependency of the output values, and
perhaps behavior, on nargout:

    [C,IB] = JOIN(...) returns an index vector IB, where JOIN constructs C by
    horizontally concatenating A(:,LEFTVARS) and B(IB,RIGHTVARS).

	...

	[C,IA,IB] = JOIN(A, B, 'Type',TYPE, ...) returns index vectors IA and IB
    indicating the correspondence between observations in C and those in A and
    B.

Whaaaat?! IS the second output an index into A, or into B? Or does the
behavior of join actually change based on nargout? Does it produce
left outer joins when nargout is 2 and inner joins when nargout is 3?

Also, you can't do the join using no key bariables. (Yes, you want to
do this often, to get the cartesian product of two sets. the
dataset equivalent of `meshgrid`.)
