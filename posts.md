---
layout: page
title: Tutorials & Notes
permalink: /posts/
lead: Practical scientific-computing workflows, theory notes behind the tools, and references for the optics lab.
description: Tutorials, scientific-computing workflows, theory notes and lab references from Indra's Lab.
---

<div class="card-grid">
{%- assign items = site.data.tools | where: "category", "learn" -%}
{%- for t in items -%}{% include tool-card.html tool=t %}{%- endfor -%}
{%- assign ref = site.data.tools | where: "id", "reference" | first -%}
{% include tool-card.html tool=ref %}
</div>

## All articles

<ul class="article-list">
{%- assign notes = site.pages | where_exp: "p", "p.path contains 'posts/'" | sort: "title" -%}
{%- for p in notes -%}
  {%- if p.title -%}
  <li><a href="{{ p.url | relative_url }}">{{ p.title }}</a>{% if p.description %}<span>{{ p.description }}</span>{% endif %}</li>
  {%- endif -%}
{%- endfor -%}
</ul>

The theory sections of the [TMM solver](../tools/tmm/#tm-theory), the [Feynman diagram generator](../tools/dsfd/) and each [calculator](../tools/) double as short notes, with equations, definitions and references.

{% include preview-scripts.html %}
