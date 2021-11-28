# Notes

## Goals

* Must have
	* Site is available
	* Modern tools
	* Modern look and feel
	* Source code available
	* Blog
	* Trusted users can update
	* Separation of content and style
		* Themes separate
* Nice to have
	* Mobile-friendly
	* Tags for papers, protocols
	* Share
	* Search
	* Analytics

## Approach

* GitHub Pages
* Jekyll + Bootstrap (and Jekyll Bootstrap)
* Sass support for CSS
* Google Fonts

## What I'm Currently Blocked On

Issues are tracked at GitHub [github-issues].

[github-issues]: https://github.com/ewallace/ewallace.github.io/issues

* Many style choices need to be updated.
* Bottom bar
* No Projects page
* Content
	* Lab members
	* Publications
		* Extremely fast and incredibly close
		* Grad school papers
	* Protocols
	

        <p>&copy; {{ site.time | date: '%Y' }} {{ site.author.name }}
          with help from <a href="http://jekyllbootstrap.com" target="_blank" title="The Definitive Jekyll Blogging Framework">Jekyll Bootstrap</a>
          and <a href="http://twitter.github.com/bootstrap/" target="_blank">Twitter Bootstrap</a>
        </p>


This site has useful stuff on social media links etc.
http://erjjones.github.io/blog/How-I-built-my-blog-in-one-day/

<a href="https://twitter.com/share" class="twitter-share-button">Tweet</a>
<script>
	!function(d,s,id){
		var js,fjs=d.getElementsByTagName(s)[0],p=/^http:/.test(d.location)?'http':'https';
		if(!d.getElementById(id)){
			js=d.createElement(s);
			js.id=id;
			js.src=p+'://platform.twitter.com/widgets.js';
			fjs.parentNode.insertBefore(js,fjs);
		}
	}(document, 'script', 'twitter-wjs');
</script>


## Layout

I'm inspired by [Trevor Bedford]'s clean, functional [site][1]. I intend to copy it---part of Trevor's intent.

* Home
	* Research
	* Papers
	* Team
	* Blog
	* Protocols
	* Links
	* Misc
		* About
		* Join
		* Contact

Each of the major areas is structured like a blog. For example:

	research/
		index.html
		_posts/
			2015-04-18-heat-shock-aggregation.md
	papers/
		index.html
		_posts/
			2014-12-09-queuosine.md
			2015-04-18-mrna-protein.md
	protocols/
		index.html
		_posts/
			2015-04-18-yeast-lysis.md
	team/
		index.html
		_posts/
			2011-10-01-d-allan-drummond.md
			2011-10-01-edward-wallace.md
	news/
		index.html
		_posts/
			2015-04-09-pg-paper-accepted.md



## To do

* Links/collaborators page
* Home page
	* Picture
	* Links to About, Join, Contact
* Style
    * Colour, layout, distinctive
* About page
* Twitter/Disqus...


## Done

* Figure out how top navbar works

## Hosting

Name registration at [NameCheap]. Custom domain via [GitHub Pages].
[NameCheap]: https://www.namecheap.com

[Trevor Bedford]: http://bedford.io/team/trevor-bedford/
[1]: http://bedford.io

