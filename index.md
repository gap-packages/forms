---
layout: default
---

# GAP Package {{site.data.package.name}}

{{site.data.package.abstract}}

The current version of this package is version {{site.data.package.version}}, released on {{site.data.package.date}}.
For more information, please refer to [the package manual]({{site.data.package.doc-html}}).
There is also a [README](README.html) file.

# An example of how it works

gap> LoadPackage("forms");<br>
true<br>
gap> gf := GF(8);<br>
GF(2^3)<br>
gap> r := PolynomialRing( gf, 3 );<br>
PolynomialRing(..., [ x_1, x_2, x_3 ])<br>
gap> poly := r.1^2 + r.2 * r.3;<br>
x_1^2+x_2\*x_3<br>
gap> form := QuadraticFormByPolynomial( poly, r );<br>
&lt; quadratic form ><br>
gap> Display( form );<br>
Quadratic form<br>
Gram Matrix:<br>
 1 . .<br>
 . . 1<br>
 . . .<br>
Polynomial: x_1^2+x_2\*x_3<br>
<br>
gap> IsDegenerateForm( form );<br>
#I  Testing degeneracy of the \*associated bilinear form\*<br>
true<br>
gap> IsSingularForm( form );<br>
false<br>
gap> WittIndex( form );<br>
1<br>
gap> IsParabolicForm( form );<br>
true<br>
gap> form;<br>
< non-singular parabolic quadratic form ><br>

## Dependencies

This package requires GAP version {{site.data.package.GAP}}
{% if site.data.package.needed-pkgs %}
The following other GAP packages are needed:
{% for pkg in site.data.package.needed-pkgs %}
- {% if pkg.url %}<a href="{{ pkg.url }}">{{ pkg.name }}</a> {% else %}{{ pkg.name }} {% endif %}
  {{- pkg.version -}}
{% endfor %}
{% endif %}
{% if site.data.package.suggested-pkgs %}
The following additional GAP packages are not required, but suggested:
{% for pkg in site.data.package.suggested-pkgs %}
- {% if pkg.url %}<a href="{{ pkg.url }}">{{ pkg.name }}</a> {% else %}{{ pkg.name }} {% endif %}
  {{- pkg.version -}}
{% endfor %}
{% endif %}


## Author{% if site.data.package.authors.size != 1 %}s{% endif %}
{% for person in site.data.package.authors %}
 {% if person.url %}<a href="{{ person.url }}">{{ person.name }}</a>{% else %}{{ person.name }}{% endif %}
 {%- if forloop.last -%}.{% else %}, {%- endif -%}
{% endfor %}

{% if site.data.package.contributors and site.data.package.contributors.size > 0 %}
## Contributor{% if site.data.package.contributors.size != 1 %}s{% endif %}
 {% for person in site.data.package.contributors %}
  {% if person.url %}<a href="{{ person.url }}">{{ person.name }}</a>{% else %}{{ person.name }}{% endif %}
  {%- if forloop.last -%}.{% else %}, {%- endif -%}
 {% endfor %}
{% endif %}

{% if site.data.package.citeas %}
## Citing

Please, cite this package as

{{site.data.package.citeas}}

You can get more info by typing `Cite("{{ site.data.package.name }}");` in the gap prompt.

{% include button-bibtex.html %}

{% endif %}


{% if site.github.issues_url %}
## Feedback

For bug reports, feature requests and suggestions, please use the
[issue tracker]({{site.github.issues_url}}).
{% endif %}
