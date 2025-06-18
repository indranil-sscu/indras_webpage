---
layout: page
title: Posts
permalink: /posts/
---

## Project Notes

Below are all the posts (subpages) from the `posts/` directory.

<ul>
  {% for page in site.pages %}
    {% if page.path contains 'posts/' and page.title %}
      <li><a href="{{ page.url | relative_url }}">{{ page.title }}</a></li>
    {% endif %}
  {% endfor %}
</ul>