I'm working with someone, and they asked for some of my intermediate
data. In the interests of what I thought would be maximum
interoperability with whatever data analysis system they preferred, I
gave them a `.csv` file.

It looked something like this:

    "target_spacing","subject","sensitivity","direction_content"
    2.0943951023932,"as",4.38700587042417,0.1
    2.0943951023932,"cj",3.21032026169591,0.1
    2.0943951023932,"jb",6.29391042968021,0.1
    2.0943951023932,"je",-0.0328585668332463,0.1
    2.0943951023932,"ko",-0.431840201907717,0.1
    2.0943951023932,"mc",-4.49124425170629,0.1
    2.0943951023932,"ml",-0.367323617983639,0.1
    2.0943951023932,"nj",2.8780230556088,0.1

Simple enough, right? For example, if you wanted to do a scatterplot
of "sensitivity" versus "target_spacing," symbols coded by subject and
colors coded by direction content, you do something like:

    library(ggplot2)
    data &lt;- read.csv('dataset.csv')
    qplot(data=data, target_spacing, sensitivity, 
          color=folded_direction_content, pch=subject, 
    	  geom="point")

On the other hand, the MATLAB script I got back from this person
looked more like

    [ num, txt, data]=xlsread('datafile.csv', 'A2:Q192');
    
    sublist= unique(txt(:, 5));
    spacinglist=unique(num(:, 4));
    colorlist=hot(length(dirlist)+3);
    symlist={'o', 's' '*' 'x' 'd' '^' 'v' '&gt;' '&lt;' 'p' 'h' 'v' '&gt;' '&lt;' 'p' 'h'};

    for s=1:length(sublist)
        for d=1:length(dirlist)
            for sp=1:length(spacinglist)
                ind=find( num(:, 3)==dirlist(d) & num(:, 4)==spacinglist(sp) & strcmp(txt(:,5),sublist{s}));
                if ~isempty(ind)
                    plot(spacinglist(sp), num(ind, 6), [symlist{s}], ...
                     'MarkerSize', 5, 'Color', colorlist(d, :),'MarkerFaceColor', ...
                     colorlist(d, :));
                     hold on
                end
            end
        end
    end

Well, this is fairly typical for MATLAB code as it is found in the wild. Give people
matrices and they use explicit numeric indices for everything. Aside
from the difficulties of making marker shape and color covary with
dimensions of the data, you have to open up the file in Excel and
count along its columns to see what variable they think they're
plotting (and after some head-scratching it turns out they weren't,
actually, plotting what they thought.)

One of the very useful features of R is that you can assign names
almost everywhere you would use an index. So, you never have to worry
about whether column 4 is "target\_spacing" or something else. You
just say "target\_spacing".

For example, let's say you have some nice rectilinear data, like this
cross-tabulation of hair color, and eye color, and sex in a group of
students:

    data &lt;- data(HairEyeColor)
    &gt; HairEyeColor
    , , Sex = Male
    
           Eye
    Hair    Brown Blue Hazel Green
      Black    32   11    10     3
      Brown    53   50    25    15
      Red      10   10     7     7
      Blond     3   30     5     8

    , , Sex = Female

           Eye
    Hair    Brown Blue Hazel Green
      Black    36    9     5     2
      Brown    66   34    29    14
      Red      16    7     7     7
      Blond     4   64     5     8

This is just a 3-D array, like Matlab's 3-D arrays (Interestingly,
Matlab only added multi-D arrays _after_ someone got fed up with the
lack of them and went off to write the Numeric package for Python.)
And as an aside, NumPy and R have consistent rules for indexing in N-D
(where N can be 1, 2, 3, or more), while MATLAB forgets about 1
dimensional arrays entirely and as for consistency,
[utterly screwed it up](/2009/08/01/for-your-thoughts-a-question-in-three-parts/).

Ahem. As I was saying, unlike an array in Matlab, arrays in R can have
nice, human-interpretable names attached to their rows, columns, and
slices. You can see them in the printout above, or get and set them
explicitly with `dimnames`:

    &gt; dimnames(HairEyeColor)
    $Hair
    [1] "Black" "Brown" "Red"   "Blond"
    
    $Eye
    [1] "Brown" "Blue"  "Hazel" "Green"
    
    $Sex
    [1] "Male"   "Female"

An array with dimnames, allows you to access elements by name, not
number. So if you want to slice just the blond, brown-eyed people
in this sample, you can just say:

    &gt; HairEyeColor['Blond', 'Brown',]
      Male Female 
         3      4 

That's the same as writing `HairEyeColor[4, 1,]`, only you can
actually see what it's trying to accomplish.

Now, I wish that you would be able to go a step further and write
`HairEyeColor[Eye='Brown',Hair='Blue',]`, and not worry about which
order the dimensions come in, but R's not perfect. Just
useful. Actually, you _can_ do that sort of thing with PANDAS, a Python
library billing itself as "[R's `data.frame` on
steroids](http://pyvideo.org/video/696/pandas-powerful-data-analysis-tools-for-python)."

Meanwhile, if you pay an additional tithe to the Mathworks, you can
get the Statistics toolbox, whose "dataset" class is more or less R's
`data.frame` with hyponatremia. (No 'NA', you see.)

Anyway, if you ever ask me to remember that "Female" is 1 (in this
dataset) and "Hazel" is 3, well, look, remembering arbitrary
correspondences between names and numbers is something humans are just
really bad at and computers are very good at, OK? If you're writing
analysis scripts and you find yourself flipping back to the speadcheet
to count columns... just don't. Why would you do a job the computer
should be doing for you?

Okay. Before I had the first gin and tonic and decided to cover a
topic or two on the syllabus of Stuff That's In Every Useful
Programming Language Except MATLAB 101, I had this script someone sent
me, that read in some data from a CSV file I'd sent them. And they
were using numeric indices into the data because they had just loaded
in the data as an array, using `xlsread`, which which doesn't do
anything useful about column headers. But you ought to be able to load
each column of data into a separate field of a struct, use the column
headers as struct field names, and refer to them by name that way, you
know, and that'd be doing pretty good for MATLAB. So I was planning to
tweak this code and send it back with a note about "here's a nice way
to do it better and let the computer take more of your mental load" 
(this person teaches a course on MATLAB for
scientists, you see, so I want to slightly reduce the fucking brain
damage that gets propagated out into the academic world.)

All you'd have to do is, instead of reading a CSV as a matrix, use the
function that reads from a CSV file and uses the column headings to
assign the fields in a struct. You know, that function. The one that
does the single bleeding obvious thing to to with CSV files. You know,
the CSV-reading function. I mean for all I rant about it, people get
work done with MATLAB. It's just impossible that Matlab can't read the
world's most ubiquitous tabular data format usefully. Right?

Well, let's try it.

The first thing I find is `csvread`. Aside from being deprecated
according to the documentation, there's another problem in that it
only reads numeric data. Now, some of the columns in my file have
things like a human observer's initials, or categorical data that's
better expressed to humans with labels like "left" or "right" rather
than trying to remember which one of those correspond to zero and 1.
(R has a built in "factor" data type to handle categorically
enumerable data, while MATLAB has bupkis.) So, `csvread` can't cut it,
because it only handles numeric data. Same problem with `dlmread`.

Next up we have `xlsread`. That's what my collaborator used to begin
with. Maybe it has an option to get the column names. Well, it won't
even read the file on my MacBook. Nor the Linux cluster we have in our
lab. Ah, see, `xlsread` only reads a CSV file if it can farm its work
out via a goddamn COM call to a motherfucking installed copy of Microsoft
Excel, and it only knows how to do that on %$)@&amp;#%%..% Windows. And, even if
my computers met those conditions, `xlsread` doesn't read a
file with more than 2^16 rows. Man, I've got more than 2^16 rows
sitting here just from asking people to look at things and press
buttons. Lord help me if I ever have a real _dataset._

CSV, you know, one of the world's most ubiquitous, plain-text,
human-readable file formats.

What next? There's `importdata` which purports to DWIM the reading of
tabular data. And there's the "Data Import Wizard" which just turns
out to be a wrapper for `importdata`.

Except `importdata` doesn't handle the quoting conventions in CSV
files. Even if that weren't a problem, it doesn't work at all. It
detects that there's a header row but it doesn't actually give me the
field names--why? Some experimentation reveals that it's, again,
completely incapable of handling non-numeric data in columns -- even
though it purports to put out separate 'data' and 'textdata' results!
Here's how 'importdata' mangles a perfectly straightforward file:

    &gt;&gt; type testfile.txt
    
    Height,Width,Depth,Weight,Label,Age,Speed
    95.01,76.21,61.54,40.57,Charlie,20.28,1.53
    23.11,45.65,79.19,93.55,Echo,19.87,74.68
    60.68,1.85,92.18,91.69,Delta,60.38,44.51
    48.60,82.14,73.82,41.03,Alpha,27.22,93.18
    89.13,44.47,17.63,89.36,Romeo,19.88,46.60
    &gt;&gt; [a delim nhreaderlines] = importdata('testfile.txt')
    a = 
            data: [5x2 double]
        textdata: {6x7 cell}
    delim =
    ,
    nheaderlines =
         1
    &gt;&gt; a.textdata
    ans = 
      Columns 1 through 6
        'Height'    'Width'    'Depth'    'Weight'    'Label'      'Age'
        '95.01'     '76.21'    '61.54'    '40.57'     'Charlie'    ''   
        '23.11'     '45.65'    '79.19'    '93.55'     'Echo'       ''   
        '60.68'     '1.85'     '92.18'    '91.69'     'Delta'      ''   
        '48.60'     '82.14'    '73.82'    '41.03'     'Alpha'      ''   
        '89.13'     '44.47'    '17.63'    '89.36'     'Romeo'      ''   
      Column 7
        'Speed'
        ''     
        ''     
        ''     
        ''     
        ''     
    &gt;&gt; a.data
    ans =
       20.2800    1.5300
       19.8700   74.6800
       60.3800   44.5100
       27.2200   93.1800
       19.8800   46.6000

So, it detects the delimiter and the single header row, but it doesn't
give back column names...why? The 'textdata' is full of perfectly
reasonable numeric strings that haven't been converted into, y'know,
numbers, but some of them have been blanked out. The 'data' pulled out
a _minority_ of the numeric data but gives you no idea _which_ columns
it pulled out for you, and that's the _best_ that I've seen so far.

The File exchange was not helpful. `csv2struct` was just a wrapper for
`xlsread` (requires Excel on Windows, limited to 65,535 rows.)
`txt2mat` claimed to be ultimately versatile and able to handle
mixed-datatype CSV files, but rejected everything I gave to it, unless
I threw enough regular expressions at its options that I might as well
have written my own CSV parser.

So I ended up writing my own fucking CSV parser.  And at several points
I got waylaid by things like:

* `textscan` will skip over blank fields when they occur at the end of
  the line (which are common in actual data), and if there isn't a
  consistent number of fields (that it doesn't skip) per line, it will
  cheerfully forget line boundaries for you. So you need to do
  conversion line-at-a-time.

* There's no good way to convert a cell array of strings to numbers,
  `str2double` tries, but it outputs NaNs whenever there's an empty
  string, or anything else it can't convert. So there's no way to tell
  whether some converted value is NaN because the file said "NaN"
  versus whether some value is NaN because the file said "booger
  police." See, the thing is, NaN is a value in the IEEE 754 system
  that is used to represent undefined values, invalid operations, and
  the like. The purpose of NaN is to signal Problems that happened
  with your arithmetic. *NaN is not a missing value marker, unles you
  really want something to to obscure places where your math is going
  wrong.* (This is why R and PANDAS allow explicitly missing -- not
  NaN -- values, in any vector -- not just floats.)

* MATLAB's version of `sscanf` can't apply an alternate locale. If you ever
  interact with people outside the US you will see some CSV type files
  with conventions like: "," for the decimal separator, "." for a
  thousands place separator, and ";" for the field delimiter. That's
  why the whole LOCALE facility in the C standard library exists. R
  provides the ability to set the locale for the purposes of reading
  such a file; whereas MATLAB's documentation explicitly forbids
  setting the locale, even in a MEX function.
  
* Speaking of MEX functions, I might have gotten this done faster if I
  had gone that route and done the parser in C/flex/bison like a grownup, instead of
  expecting MATLAB to be any help at all in doing stuff like 
  converting several strings to numbers.

So, as you see, reading a CSV file into MATLAB entails a whole lot of
bullshit.

By comparison, here's how you read that exact same file in R.

    &gt; x = read.csv("testfile.txt")
    &gt; x
      Height Width Depth Weight   Label   Age Speed
    1  95.01 76.21 61.54  40.57 Charlie 20.28  1.53
    2  23.11 45.65 79.19  93.55    Echo 19.87 74.68
    3  60.68  1.85 92.18  91.69   Delta 60.38 44.51
    4  48.60 82.14 73.82  41.03   Alpha 27.22 93.18
    5  89.13 44.47 17.63  89.36   Romeo 19.88 46.60
    &gt; class(x$Width)
    [1] "numeric"
    &gt; class(x$Name)
    [1] "factor"

You see? Ask R to read in a table, and it makes a good guess at the
appropriate data types and headers, and you can refer to the
components by their actual names. This stuff just _isn't so hard._
