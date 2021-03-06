




<!DOCTYPE html>
<html class="   ">
  <head prefix="og: http://ogp.me/ns# fb: http://ogp.me/ns/fb# object: http://ogp.me/ns/object# article: http://ogp.me/ns/article# profile: http://ogp.me/ns/profile#">
    <meta charset='utf-8'>
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    
    
    <title>cmssw/GeneratorInterface/RivetInterface/src/CMS_SMP_12_017.cc at 53XZjets · fabiocos/cmssw · GitHub</title>
    <link rel="search" type="application/opensearchdescription+xml" href="/opensearch.xml" title="GitHub" />
    <link rel="fluid-icon" href="https://github.com/fluidicon.png" title="GitHub" />
    <link rel="apple-touch-icon" sizes="57x57" href="/apple-touch-icon-114.png" />
    <link rel="apple-touch-icon" sizes="114x114" href="/apple-touch-icon-114.png" />
    <link rel="apple-touch-icon" sizes="72x72" href="/apple-touch-icon-144.png" />
    <link rel="apple-touch-icon" sizes="144x144" href="/apple-touch-icon-144.png" />
    <meta property="fb:app_id" content="1401488693436528"/>

      <meta content="@github" name="twitter:site" /><meta content="summary" name="twitter:card" /><meta content="fabiocos/cmssw" name="twitter:title" /><meta content="cmssw - CMS Offline Software" name="twitter:description" /><meta content="https://avatars2.githubusercontent.com/u/4058194?s=400" name="twitter:image:src" />
<meta content="GitHub" property="og:site_name" /><meta content="object" property="og:type" /><meta content="https://avatars2.githubusercontent.com/u/4058194?s=400" property="og:image" /><meta content="fabiocos/cmssw" property="og:title" /><meta content="https://github.com/fabiocos/cmssw" property="og:url" /><meta content="cmssw - CMS Offline Software" property="og:description" />

    <link rel="assets" href="https://assets-cdn.github.com/">
    <link rel="conduit-xhr" href="https://ghconduit.com:25035">
    

    <meta name="msapplication-TileImage" content="/windows-tile.png" />
    <meta name="msapplication-TileColor" content="#ffffff" />
    <meta name="selected-link" value="repo_source" data-pjax-transient />
      <meta name="google-analytics" content="UA-3769691-2">

    <meta content="collector.githubapp.com" name="octolytics-host" /><meta content="collector-cdn.github.com" name="octolytics-script-host" /><meta content="github" name="octolytics-app-id" /><meta content="808DF556:01FE:9027CCD:53A07B01" name="octolytics-dimension-request_id" />
    

    
    
    <link rel="icon" type="image/x-icon" href="https://assets-cdn.github.com/favicon.ico" />


    <meta content="authenticity_token" name="csrf-param" />
<meta content="omEY8noD9jUHEqh8P0q2HABcWFdj/Fa+X0E/uY21Kxo/jqgi7TiBndSoCghViRhz68M0bhk3YyFHJY02QauG8g==" name="csrf-token" />

    <link href="https://assets-cdn.github.com/assets/github-bcdb9eb20e01413e751df455ef4c7bf63158bec4.css" media="all" rel="stylesheet" type="text/css" />
    <link href="https://assets-cdn.github.com/assets/github2-bc74bff9146dfddbac38e3061835e070dedb4d3b.css" media="all" rel="stylesheet" type="text/css" />
    


    <meta http-equiv="x-pjax-version" content="3ff95f6c823e42a9e1ea5fec8b73ae56">

      
  <meta name="description" content="cmssw - CMS Offline Software" />

  <meta content="4058194" name="octolytics-dimension-user_id" /><meta content="fabiocos" name="octolytics-dimension-user_login" /><meta content="11260036" name="octolytics-dimension-repository_id" /><meta content="fabiocos/cmssw" name="octolytics-dimension-repository_nwo" /><meta content="true" name="octolytics-dimension-repository_public" /><meta content="true" name="octolytics-dimension-repository_is_fork" /><meta content="10969551" name="octolytics-dimension-repository_parent_id" /><meta content="cms-sw/cmssw" name="octolytics-dimension-repository_parent_nwo" /><meta content="10969551" name="octolytics-dimension-repository_network_root_id" /><meta content="cms-sw/cmssw" name="octolytics-dimension-repository_network_root_nwo" />
  <link href="https://github.com/fabiocos/cmssw/commits/53XZjets.atom" rel="alternate" title="Recent Commits to cmssw:53XZjets" type="application/atom+xml" />

  </head>


  <body class="logged_out  env-production linux vis-public fork page-blob">
    <a href="#start-of-content" tabindex="1" class="accessibility-aid js-skip-to-content">Skip to content</a>
    <div class="wrapper">
      
      
      
      


      
      <div class="header header-logged-out">
  <div class="container clearfix">

    <a class="header-logo-wordmark" href="https://github.com/">
      <span class="mega-octicon octicon-logo-github"></span>
    </a>

    <div class="header-actions">
        <a class="button primary" href="/join">Sign up</a>
      <a class="button signin" href="/login?return_to=%2Ffabiocos%2Fcmssw%2Fblob%2F53XZjets%2FGeneratorInterface%2FRivetInterface%2Fsrc%2FCMS_SMP_12_017.cc">Sign in</a>
    </div>

    <div class="command-bar js-command-bar  in-repository">

      <ul class="top-nav">
          <li class="explore"><a href="/explore">Explore</a></li>
        <li class="features"><a href="/features">Features</a></li>
          <li class="enterprise"><a href="https://enterprise.github.com/">Enterprise</a></li>
          <li class="blog"><a href="/blog">Blog</a></li>
      </ul>
        <form accept-charset="UTF-8" action="/search" class="command-bar-form" id="top_search_form" method="get">

<div class="commandbar">
  <span class="message"></span>
  <input type="text" data-hotkey="s" name="q" id="js-command-bar-field" placeholder="Search or type a command" tabindex="1" autocapitalize="off"
    
    
      data-repo="fabiocos/cmssw"
      data-branch="53XZjets"
      data-sha="04079e6106bf05fe0762d4d6404bb5434e0905c9"
  >
  <div class="display hidden"></div>
</div>

    <input type="hidden" name="nwo" value="fabiocos/cmssw" />

    <div class="select-menu js-menu-container js-select-menu search-context-select-menu">
      <span class="minibutton select-menu-button js-menu-target" role="button" aria-haspopup="true">
        <span class="js-select-button">This repository</span>
      </span>

      <div class="select-menu-modal-holder js-menu-content js-navigation-container" aria-hidden="true">
        <div class="select-menu-modal">

          <div class="select-menu-item js-navigation-item js-this-repository-navigation-item selected">
            <span class="select-menu-item-icon octicon octicon-check"></span>
            <input type="radio" class="js-search-this-repository" name="search_target" value="repository" checked="checked" />
            <div class="select-menu-item-text js-select-button-text">This repository</div>
          </div> <!-- /.select-menu-item -->

          <div class="select-menu-item js-navigation-item js-all-repositories-navigation-item">
            <span class="select-menu-item-icon octicon octicon-check"></span>
            <input type="radio" name="search_target" value="global" />
            <div class="select-menu-item-text js-select-button-text">All repositories</div>
          </div> <!-- /.select-menu-item -->

        </div>
      </div>
    </div>

  <span class="help tooltipped tooltipped-s" aria-label="Show command bar help">
    <span class="octicon octicon-question"></span>
  </span>


  <input type="hidden" name="ref" value="cmdform">

</form>
    </div>

  </div>
</div>



      <div id="start-of-content" class="accessibility-aid"></div>
          <div class="site" itemscope itemtype="http://schema.org/WebPage">
    <div id="js-flash-container">
      
    </div>
    <div class="pagehead repohead instapaper_ignore readability-menu">
      <div class="container">
        

<ul class="pagehead-actions">


  <li>
      <a href="/login?return_to=%2Ffabiocos%2Fcmssw"
    class="minibutton with-count star-button tooltipped tooltipped-n"
    aria-label="You must be signed in to star a repository" rel="nofollow">
    <span class="octicon octicon-star"></span>
    Star
  </a>

    <a class="social-count js-social-count" href="/fabiocos/cmssw/stargazers">
      0
    </a>

  </li>

    <li>
      <a href="/login?return_to=%2Ffabiocos%2Fcmssw"
        class="minibutton with-count js-toggler-target fork-button tooltipped tooltipped-n"
        aria-label="You must be signed in to fork a repository" rel="nofollow">
        <span class="octicon octicon-repo-forked"></span>
        Fork
      </a>
      <a href="/fabiocos/cmssw/network" class="social-count">
        803
      </a>
    </li>
</ul>

        <h1 itemscope itemtype="http://data-vocabulary.org/Breadcrumb" class="entry-title public">
          <span class="repo-label"><span>public</span></span>
          <span class="mega-octicon octicon-repo"></span>
          <span class="author"><a href="/fabiocos" class="url fn" itemprop="url" rel="author"><span itemprop="title">fabiocos</span></a></span><!--
       --><span class="path-divider">/</span><!--
       --><strong><a href="/fabiocos/cmssw" class="js-current-repository js-repo-home-link">cmssw</a></strong>

          <span class="page-context-loader">
            <img alt="" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
          </span>

            <span class="fork-flag">
              <span class="text">forked from <a href="/cms-sw/cmssw">cms-sw/cmssw</a></span>
            </span>
        </h1>
      </div><!-- /.container -->
    </div><!-- /.repohead -->

    <div class="container">
      <div class="repository-with-sidebar repo-container new-discussion-timeline js-new-discussion-timeline  ">
        <div class="repository-sidebar clearfix">
            

<div class="sunken-menu vertical-right repo-nav js-repo-nav js-repository-container-pjax js-octicon-loaders">
  <div class="sunken-menu-contents">
    <ul class="sunken-menu-group">
      <li class="tooltipped tooltipped-w" aria-label="Code">
        <a href="/fabiocos/cmssw/tree/53XZjets" aria-label="Code" class="selected js-selected-navigation-item sunken-menu-item" data-hotkey="g c" data-pjax="true" data-selected-links="repo_source repo_downloads repo_commits repo_releases repo_tags repo_branches /fabiocos/cmssw/tree/53XZjets">
          <span class="octicon octicon-code"></span> <span class="full-word">Code</span>
          <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>


      <li class="tooltipped tooltipped-w" aria-label="Pull Requests">
        <a href="/fabiocos/cmssw/pulls" aria-label="Pull Requests" class="js-selected-navigation-item sunken-menu-item js-disable-pjax" data-hotkey="g p" data-selected-links="repo_pulls /fabiocos/cmssw/pulls">
            <span class="octicon octicon-git-pull-request"></span> <span class="full-word">Pull Requests</span>
            <span class='counter'>0</span>
            <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>


    </ul>
    <div class="sunken-menu-separator"></div>
    <ul class="sunken-menu-group">

      <li class="tooltipped tooltipped-w" aria-label="Pulse">
        <a href="/fabiocos/cmssw/pulse" aria-label="Pulse" class="js-selected-navigation-item sunken-menu-item" data-pjax="true" data-selected-links="pulse /fabiocos/cmssw/pulse">
          <span class="octicon octicon-pulse"></span> <span class="full-word">Pulse</span>
          <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>

      <li class="tooltipped tooltipped-w" aria-label="Graphs">
        <a href="/fabiocos/cmssw/graphs" aria-label="Graphs" class="js-selected-navigation-item sunken-menu-item" data-pjax="true" data-selected-links="repo_graphs repo_contributors /fabiocos/cmssw/graphs">
          <span class="octicon octicon-graph"></span> <span class="full-word">Graphs</span>
          <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>

      <li class="tooltipped tooltipped-w" aria-label="Network">
        <a href="/fabiocos/cmssw/network" aria-label="Network" class="js-selected-navigation-item sunken-menu-item js-disable-pjax" data-selected-links="repo_network /fabiocos/cmssw/network">
          <span class="octicon octicon-repo-forked"></span> <span class="full-word">Network</span>
          <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>
    </ul>


  </div>
</div>

              <div class="only-with-full-nav">
                

  

<div class="clone-url open"
  data-protocol-type="http"
  data-url="/users/set_protocol?protocol_selector=http&amp;protocol_type=clone">
  <h3><strong>HTTPS</strong> clone URL</h3>
  <div class="clone-url-box">
    <input type="text" class="clone js-url-field"
           value="https://github.com/fabiocos/cmssw.git" readonly="readonly">
    <span class="url-box-clippy">
    <button aria-label="copy to clipboard" class="js-zeroclipboard minibutton zeroclipboard-button" data-clipboard-text="https://github.com/fabiocos/cmssw.git" data-copied-hint="copied!" type="button"><span class="octicon octicon-clippy"></span></button>
    </span>
  </div>
</div>

  

<div class="clone-url "
  data-protocol-type="subversion"
  data-url="/users/set_protocol?protocol_selector=subversion&amp;protocol_type=clone">
  <h3><strong>Subversion</strong> checkout URL</h3>
  <div class="clone-url-box">
    <input type="text" class="clone js-url-field"
           value="https://github.com/fabiocos/cmssw" readonly="readonly">
    <span class="url-box-clippy">
    <button aria-label="copy to clipboard" class="js-zeroclipboard minibutton zeroclipboard-button" data-clipboard-text="https://github.com/fabiocos/cmssw" data-copied-hint="copied!" type="button"><span class="octicon octicon-clippy"></span></button>
    </span>
  </div>
</div>


<p class="clone-options">You can clone with
      <a href="#" class="js-clone-selector" data-protocol="http">HTTPS</a>
      or <a href="#" class="js-clone-selector" data-protocol="subversion">Subversion</a>.
  <a href="https://help.github.com/articles/which-remote-url-should-i-use" class="help tooltipped tooltipped-n" aria-label="Get help on which URL is right for you.">
    <span class="octicon octicon-question"></span>
  </a>
</p>



                <a href="/fabiocos/cmssw/archive/53XZjets.zip"
                   class="minibutton sidebar-button"
                   aria-label="Download fabiocos/cmssw as a zip file"
                   title="Download fabiocos/cmssw as a zip file"
                   rel="nofollow">
                  <span class="octicon octicon-cloud-download"></span>
                  Download ZIP
                </a>
              </div>
        </div><!-- /.repository-sidebar -->

        <div id="js-repo-pjax-container" class="repository-content context-loader-container" data-pjax-container>
          


<a href="/fabiocos/cmssw/blob/fa74cc6c537d1605638fa40f06b0409fbdc9ccaa/GeneratorInterface/RivetInterface/src/CMS_SMP_12_017.cc" class="hidden js-permalink-shortcut" data-hotkey="y">Permalink</a>

<!-- blob contrib key: blob_contributors:v21:32238f175fae39925a211cd78a20d51b -->

<p title="This is a placeholder element" class="js-history-link-replace hidden"></p>

<div class="file-navigation">
  

<div class="select-menu js-menu-container js-select-menu" >
  <span class="minibutton select-menu-button js-menu-target css-truncate" data-hotkey="w"
    data-master-branch="CMSSW_7_0_X"
    data-ref="53XZjets"
    title="53XZjets"
    role="button" aria-label="Switch branches or tags" tabindex="0" aria-haspopup="true">
    <span class="octicon octicon-git-branch"></span>
    <i>branch:</i>
    <span class="js-select-button css-truncate-target">53XZjets</span>
  </span>

  <div class="select-menu-modal-holder js-menu-content js-navigation-container" data-pjax aria-hidden="true">

    <div class="select-menu-modal">
      <div class="select-menu-header">
        <span class="select-menu-title">Switch branches/tags</span>
        <span class="octicon octicon-x js-menu-close"></span>
      </div> <!-- /.select-menu-header -->

      <div class="select-menu-filters">
        <div class="select-menu-text-filter">
          <input type="text" aria-label="Filter branches/tags" id="context-commitish-filter-field" class="js-filterable-field js-navigation-enable" placeholder="Filter branches/tags">
        </div>
        <div class="select-menu-tabs">
          <ul>
            <li class="select-menu-tab">
              <a href="#" data-tab-filter="branches" class="js-select-menu-tab">Branches</a>
            </li>
            <li class="select-menu-tab">
              <a href="#" data-tab-filter="tags" class="js-select-menu-tab">Tags</a>
            </li>
          </ul>
        </div><!-- /.select-menu-tabs -->
      </div><!-- /.select-menu-filters -->

      <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket" data-tab-filter="branches">

        <div data-filterable-for="context-commitish-filter-field" data-filterable-type="substring">


            <div class="select-menu-item js-navigation-item selected">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/fabiocos/cmssw/blob/53XZjets/GeneratorInterface/RivetInterface/src/CMS_SMP_12_017.cc"
                 data-name="53XZjets"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="53XZjets">53XZjets</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/fabiocos/cmssw/blob/70X-LHAPDFFrontier/GeneratorInterface/RivetInterface/src/CMS_SMP_12_017.cc"
                 data-name="70X-LHAPDFFrontier"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="70X-LHAPDFFrontier">70X-LHAPDFFrontier</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/fabiocos/cmssw/blob/442Zjets/GeneratorInterface/RivetInterface/src/CMS_SMP_12_017.cc"
                 data-name="442Zjets"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="442Zjets">442Zjets</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/fabiocos/cmssw/blob/CMSSW_5_3_X/GeneratorInterface/RivetInterface/src/CMS_SMP_12_017.cc"
                 data-name="CMSSW_5_3_X"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="CMSSW_5_3_X">CMSSW_5_3_X</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/fabiocos/cmssw/blob/CMSSW_6_2_X/GeneratorInterface/RivetInterface/src/CMS_SMP_12_017.cc"
                 data-name="CMSSW_6_2_X"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="CMSSW_6_2_X">CMSSW_6_2_X</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/fabiocos/cmssw/blob/CMSSW_7_0_X/GeneratorInterface/RivetInterface/src/CMS_SMP_12_017.cc"
                 data-name="CMSSW_7_0_X"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="CMSSW_7_0_X">CMSSW_7_0_X</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/fabiocos/cmssw/blob/gh-pages/GeneratorInterface/RivetInterface/src/CMS_SMP_12_017.cc"
                 data-name="gh-pages"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="gh-pages">gh-pages</a>
            </div> <!-- /.select-menu-item -->
        </div>

          <div class="select-menu-no-results">Nothing to show</div>
      </div> <!-- /.select-menu-list -->

      <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket" data-tab-filter="tags">
        <div data-filterable-for="context-commitish-filter-field" data-filterable-type="substring">


            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/fabiocos/cmssw/tree/CMSSW_7_0_X_2013-07-08-1400/GeneratorInterface/RivetInterface/src/CMS_SMP_12_017.cc"
                 data-name="CMSSW_7_0_X_2013-07-08-1400"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="CMSSW_7_0_X_2013-07-08-1400">CMSSW_7_0_X_2013-07-08-1400</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/fabiocos/cmssw/tree/CMSSW_7_0_X_2013-07-08-0200/GeneratorInterface/RivetInterface/src/CMS_SMP_12_017.cc"
                 data-name="CMSSW_7_0_X_2013-07-08-0200"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="CMSSW_7_0_X_2013-07-08-0200">CMSSW_7_0_X_2013-07-08-0200</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/fabiocos/cmssw/tree/CMSSW_7_0_X_2013-07-07-1400/GeneratorInterface/RivetInterface/src/CMS_SMP_12_017.cc"
                 data-name="CMSSW_7_0_X_2013-07-07-1400"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="CMSSW_7_0_X_2013-07-07-1400">CMSSW_7_0_X_2013-07-07-1400</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/fabiocos/cmssw/tree/CMSSW_7_0_X_2013-07-07-0200/GeneratorInterface/RivetInterface/src/CMS_SMP_12_017.cc"
                 data-name="CMSSW_7_0_X_2013-07-07-0200"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="CMSSW_7_0_X_2013-07-07-0200">CMSSW_7_0_X_2013-07-07-0200</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/fabiocos/cmssw/tree/CMSSW_7_0_X_2013-07-06-1400/GeneratorInterface/RivetInterface/src/CMS_SMP_12_017.cc"
                 data-name="CMSSW_7_0_X_2013-07-06-1400"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="CMSSW_7_0_X_2013-07-06-1400">CMSSW_7_0_X_2013-07-06-1400</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/fabiocos/cmssw/tree/CMSSW_7_0_X_2013-07-06-0200/GeneratorInterface/RivetInterface/src/CMS_SMP_12_017.cc"
                 data-name="CMSSW_7_0_X_2013-07-06-0200"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="CMSSW_7_0_X_2013-07-06-0200">CMSSW_7_0_X_2013-07-06-0200</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/fabiocos/cmssw/tree/CMSSW_7_0_X_2013-07-05-1400/GeneratorInterface/RivetInterface/src/CMS_SMP_12_017.cc"
                 data-name="CMSSW_7_0_X_2013-07-05-1400"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="CMSSW_7_0_X_2013-07-05-1400">CMSSW_7_0_X_2013-07-05-1400</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/fabiocos/cmssw/tree/CMSSW_7_0_X_2013-07-05-0200/GeneratorInterface/RivetInterface/src/CMS_SMP_12_017.cc"
                 data-name="CMSSW_7_0_X_2013-07-05-0200"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="CMSSW_7_0_X_2013-07-05-0200">CMSSW_7_0_X_2013-07-05-0200</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/fabiocos/cmssw/tree/CMSSW_7_0_X_2013-07-04-1400/GeneratorInterface/RivetInterface/src/CMS_SMP_12_017.cc"
                 data-name="CMSSW_7_0_X_2013-07-04-1400"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="CMSSW_7_0_X_2013-07-04-1400">CMSSW_7_0_X_2013-07-04-1400</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/fabiocos/cmssw/tree/CMSSW_7_0_X_2013-07-04-0200/GeneratorInterface/RivetInterface/src/CMS_SMP_12_017.cc"
                 data-name="CMSSW_7_0_X_2013-07-04-0200"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="CMSSW_7_0_X_2013-07-04-0200">CMSSW_7_0_X_2013-07-04-0200</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/fabiocos/cmssw/tree/CMSSW_7_0_X_2013-07-03-1400/GeneratorInterface/RivetInterface/src/CMS_SMP_12_017.cc"
                 data-name="CMSSW_7_0_X_2013-07-03-1400"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="CMSSW_7_0_X_2013-07-03-1400">CMSSW_7_0_X_2013-07-03-1400</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/fabiocos/cmssw/tree/CMSSW_7_0_X_2013-07-03-0200/GeneratorInterface/RivetInterface/src/CMS_SMP_12_017.cc"
                 data-name="CMSSW_7_0_X_2013-07-03-0200"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="CMSSW_7_0_X_2013-07-03-0200">CMSSW_7_0_X_2013-07-03-0200</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/fabiocos/cmssw/tree/CMSSW_7_0_X_2013-07-02-1400/GeneratorInterface/RivetInterface/src/CMS_SMP_12_017.cc"
                 data-name="CMSSW_7_0_X_2013-07-02-1400"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="CMSSW_7_0_X_2013-07-02-1400">CMSSW_7_0_X_2013-07-02-1400</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/fabiocos/cmssw/tree/CMSSW_7_0_X_2013-07-02-1100/GeneratorInterface/RivetInterface/src/CMS_SMP_12_017.cc"
                 data-name="CMSSW_7_0_X_2013-07-02-1100"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="CMSSW_7_0_X_2013-07-02-1100">CMSSW_7_0_X_2013-07-02-1100</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/fabiocos/cmssw/tree/CMSSW_7_0_X_2013-07-02-0200/GeneratorInterface/RivetInterface/src/CMS_SMP_12_017.cc"
                 data-name="CMSSW_7_0_X_2013-07-02-0200"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="CMSSW_7_0_X_2013-07-02-0200">CMSSW_7_0_X_2013-07-02-0200</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/fabiocos/cmssw/tree/CMSSW_7_0_X_2013-06-26-2200/GeneratorInterface/RivetInterface/src/CMS_SMP_12_017.cc"
                 data-name="CMSSW_7_0_X_2013-06-26-2200"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="CMSSW_7_0_X_2013-06-26-2200">CMSSW_7_0_X_2013-06-26-2200</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/fabiocos/cmssw/tree/CMSSW_7_0_0_pre0/GeneratorInterface/RivetInterface/src/CMS_SMP_12_017.cc"
                 data-name="CMSSW_7_0_0_pre0"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="CMSSW_7_0_0_pre0">CMSSW_7_0_0_pre0</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/fabiocos/cmssw/tree/CMSSW_6_2_0_pre8/GeneratorInterface/RivetInterface/src/CMS_SMP_12_017.cc"
                 data-name="CMSSW_6_2_0_pre8"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="CMSSW_6_2_0_pre8">CMSSW_6_2_0_pre8</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/fabiocos/cmssw/tree/70X-201309061349/GeneratorInterface/RivetInterface/src/CMS_SMP_12_017.cc"
                 data-name="70X-201309061349"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="70X-201309061349">70X-201309061349</a>
            </div> <!-- /.select-menu-item -->
            <div class="select-menu-item js-navigation-item ">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/fabiocos/cmssw/tree/53X-sherpa2beta2-201308051537/GeneratorInterface/RivetInterface/src/CMS_SMP_12_017.cc"
                 data-name="53X-sherpa2beta2-201308051537"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="53X-sherpa2beta2-201308051537">53X-sherpa2beta2-201308051537</a>
            </div> <!-- /.select-menu-item -->
        </div>

        <div class="select-menu-no-results">Nothing to show</div>
      </div> <!-- /.select-menu-list -->

    </div> <!-- /.select-menu-modal -->
  </div> <!-- /.select-menu-modal-holder -->
</div> <!-- /.select-menu -->

  <div class="button-group right">
    <a href="/fabiocos/cmssw/find/53XZjets"
          class="js-show-file-finder minibutton empty-icon tooltipped tooltipped-s"
          data-pjax
          data-hotkey="t"
          aria-label="Quickly jump between files">
      <span class="octicon octicon-list-unordered"></span>
    </a>
    <button class="js-zeroclipboard minibutton zeroclipboard-button"
          data-clipboard-text="GeneratorInterface/RivetInterface/src/CMS_SMP_12_017.cc"
          aria-label="Copy to clipboard"
          data-copied-hint="Copied!">
      <span class="octicon octicon-clippy"></span>
    </button>
  </div>

  <div class="breadcrumb">
    <span class='repo-root js-repo-root'><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/fabiocos/cmssw/tree/53XZjets" data-branch="53XZjets" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">cmssw</span></a></span></span><span class="separator"> / </span><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/fabiocos/cmssw/tree/53XZjets/GeneratorInterface" data-branch="53XZjets" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">GeneratorInterface</span></a></span><span class="separator"> / </span><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/fabiocos/cmssw/tree/53XZjets/GeneratorInterface/RivetInterface" data-branch="53XZjets" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">RivetInterface</span></a></span><span class="separator"> / </span><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/fabiocos/cmssw/tree/53XZjets/GeneratorInterface/RivetInterface/src" data-branch="53XZjets" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">src</span></a></span><span class="separator"> / </span><strong class="final-path">CMS_SMP_12_017.cc</strong>
  </div>
</div>


  <div class="commit file-history-tease">
      <img alt="Fabio Cossutti" class="main-avatar js-avatar" data-user="4058194" height="24" src="https://avatars1.githubusercontent.com/u/4058194?s=140" width="24" />
      <span class="author"><a href="/fabiocos" rel="author">fabiocos</a></span>
      <time datetime="2014-03-07T17:51:57+01:00" is="relative-time">March 07, 2014</time>
      <div class="commit-title">
          <a href="/fabiocos/cmssw/commit/675a19b4fbc16f12e6c790fb222f7a78e57057de" class="message" data-pjax="true" title="20140306 update of CMS_SMP_12_017 data file">20140306 update of CMS_SMP_12_017 data file</a>
      </div>

    <div class="participation">
      <p class="quickstat"><a href="#blob_contributors_box" rel="facebox"><strong>1</strong>  contributor</a></p>
      
    </div>
    <div id="blob_contributors_box" style="display:none">
      <h2 class="facebox-header">Users who have contributed to this file</h2>
      <ul class="facebox-user-list">
          <li class="facebox-user-list-item">
            <img alt="Fabio Cossutti" class=" js-avatar" data-user="4058194" height="24" src="https://avatars1.githubusercontent.com/u/4058194?s=140" width="24" />
            <a href="/fabiocos">fabiocos</a>
          </li>
      </ul>
    </div>
  </div>

<div class="file-box">
  <div class="file">
    <div class="meta clearfix">
      <div class="info file-name">
        <span class="icon"><b class="octicon octicon-file-text"></b></span>
        <span class="mode" title="File Mode">file</span>
        <span class="meta-divider"></span>
          <span>633 lines (536 sloc)</span>
          <span class="meta-divider"></span>
        <span>22.973 kb</span>
      </div>
      <div class="actions">
        <div class="button-group">
              <a class="minibutton disabled tooltipped tooltipped-w" href="#"
                 aria-label="You must be signed in to make or propose changes">Edit</a>
          <a href="/fabiocos/cmssw/raw/53XZjets/GeneratorInterface/RivetInterface/src/CMS_SMP_12_017.cc" class="button minibutton " id="raw-url">Raw</a>
            <a href="/fabiocos/cmssw/blame/53XZjets/GeneratorInterface/RivetInterface/src/CMS_SMP_12_017.cc" class="button minibutton js-update-url-with-hash">Blame</a>
          <a href="/fabiocos/cmssw/commits/53XZjets/GeneratorInterface/RivetInterface/src/CMS_SMP_12_017.cc" class="button minibutton " rel="nofollow">History</a>
        </div><!-- /.button-group -->
          <a class="minibutton danger disabled empty-icon tooltipped tooltipped-w" href="#"
             aria-label="You must be signed in to make or propose changes">
          Delete
        </a>
      </div><!-- /.actions -->
    </div>
      
  <div class="blob-wrapper data type-c js-blob-data">
       <table class="file-code file-diff tab-size-8">
         <tr class="file-code-line">
           <td class="blob-line-nums">
             <span id="L1" rel="#L1">1</span>
<span id="L2" rel="#L2">2</span>
<span id="L3" rel="#L3">3</span>
<span id="L4" rel="#L4">4</span>
<span id="L5" rel="#L5">5</span>
<span id="L6" rel="#L6">6</span>
<span id="L7" rel="#L7">7</span>
<span id="L8" rel="#L8">8</span>
<span id="L9" rel="#L9">9</span>
<span id="L10" rel="#L10">10</span>
<span id="L11" rel="#L11">11</span>
<span id="L12" rel="#L12">12</span>
<span id="L13" rel="#L13">13</span>
<span id="L14" rel="#L14">14</span>
<span id="L15" rel="#L15">15</span>
<span id="L16" rel="#L16">16</span>
<span id="L17" rel="#L17">17</span>
<span id="L18" rel="#L18">18</span>
<span id="L19" rel="#L19">19</span>
<span id="L20" rel="#L20">20</span>
<span id="L21" rel="#L21">21</span>
<span id="L22" rel="#L22">22</span>
<span id="L23" rel="#L23">23</span>
<span id="L24" rel="#L24">24</span>
<span id="L25" rel="#L25">25</span>
<span id="L26" rel="#L26">26</span>
<span id="L27" rel="#L27">27</span>
<span id="L28" rel="#L28">28</span>
<span id="L29" rel="#L29">29</span>
<span id="L30" rel="#L30">30</span>
<span id="L31" rel="#L31">31</span>
<span id="L32" rel="#L32">32</span>
<span id="L33" rel="#L33">33</span>
<span id="L34" rel="#L34">34</span>
<span id="L35" rel="#L35">35</span>
<span id="L36" rel="#L36">36</span>
<span id="L37" rel="#L37">37</span>
<span id="L38" rel="#L38">38</span>
<span id="L39" rel="#L39">39</span>
<span id="L40" rel="#L40">40</span>
<span id="L41" rel="#L41">41</span>
<span id="L42" rel="#L42">42</span>
<span id="L43" rel="#L43">43</span>
<span id="L44" rel="#L44">44</span>
<span id="L45" rel="#L45">45</span>
<span id="L46" rel="#L46">46</span>
<span id="L47" rel="#L47">47</span>
<span id="L48" rel="#L48">48</span>
<span id="L49" rel="#L49">49</span>
<span id="L50" rel="#L50">50</span>
<span id="L51" rel="#L51">51</span>
<span id="L52" rel="#L52">52</span>
<span id="L53" rel="#L53">53</span>
<span id="L54" rel="#L54">54</span>
<span id="L55" rel="#L55">55</span>
<span id="L56" rel="#L56">56</span>
<span id="L57" rel="#L57">57</span>
<span id="L58" rel="#L58">58</span>
<span id="L59" rel="#L59">59</span>
<span id="L60" rel="#L60">60</span>
<span id="L61" rel="#L61">61</span>
<span id="L62" rel="#L62">62</span>
<span id="L63" rel="#L63">63</span>
<span id="L64" rel="#L64">64</span>
<span id="L65" rel="#L65">65</span>
<span id="L66" rel="#L66">66</span>
<span id="L67" rel="#L67">67</span>
<span id="L68" rel="#L68">68</span>
<span id="L69" rel="#L69">69</span>
<span id="L70" rel="#L70">70</span>
<span id="L71" rel="#L71">71</span>
<span id="L72" rel="#L72">72</span>
<span id="L73" rel="#L73">73</span>
<span id="L74" rel="#L74">74</span>
<span id="L75" rel="#L75">75</span>
<span id="L76" rel="#L76">76</span>
<span id="L77" rel="#L77">77</span>
<span id="L78" rel="#L78">78</span>
<span id="L79" rel="#L79">79</span>
<span id="L80" rel="#L80">80</span>
<span id="L81" rel="#L81">81</span>
<span id="L82" rel="#L82">82</span>
<span id="L83" rel="#L83">83</span>
<span id="L84" rel="#L84">84</span>
<span id="L85" rel="#L85">85</span>
<span id="L86" rel="#L86">86</span>
<span id="L87" rel="#L87">87</span>
<span id="L88" rel="#L88">88</span>
<span id="L89" rel="#L89">89</span>
<span id="L90" rel="#L90">90</span>
<span id="L91" rel="#L91">91</span>
<span id="L92" rel="#L92">92</span>
<span id="L93" rel="#L93">93</span>
<span id="L94" rel="#L94">94</span>
<span id="L95" rel="#L95">95</span>
<span id="L96" rel="#L96">96</span>
<span id="L97" rel="#L97">97</span>
<span id="L98" rel="#L98">98</span>
<span id="L99" rel="#L99">99</span>
<span id="L100" rel="#L100">100</span>
<span id="L101" rel="#L101">101</span>
<span id="L102" rel="#L102">102</span>
<span id="L103" rel="#L103">103</span>
<span id="L104" rel="#L104">104</span>
<span id="L105" rel="#L105">105</span>
<span id="L106" rel="#L106">106</span>
<span id="L107" rel="#L107">107</span>
<span id="L108" rel="#L108">108</span>
<span id="L109" rel="#L109">109</span>
<span id="L110" rel="#L110">110</span>
<span id="L111" rel="#L111">111</span>
<span id="L112" rel="#L112">112</span>
<span id="L113" rel="#L113">113</span>
<span id="L114" rel="#L114">114</span>
<span id="L115" rel="#L115">115</span>
<span id="L116" rel="#L116">116</span>
<span id="L117" rel="#L117">117</span>
<span id="L118" rel="#L118">118</span>
<span id="L119" rel="#L119">119</span>
<span id="L120" rel="#L120">120</span>
<span id="L121" rel="#L121">121</span>
<span id="L122" rel="#L122">122</span>
<span id="L123" rel="#L123">123</span>
<span id="L124" rel="#L124">124</span>
<span id="L125" rel="#L125">125</span>
<span id="L126" rel="#L126">126</span>
<span id="L127" rel="#L127">127</span>
<span id="L128" rel="#L128">128</span>
<span id="L129" rel="#L129">129</span>
<span id="L130" rel="#L130">130</span>
<span id="L131" rel="#L131">131</span>
<span id="L132" rel="#L132">132</span>
<span id="L133" rel="#L133">133</span>
<span id="L134" rel="#L134">134</span>
<span id="L135" rel="#L135">135</span>
<span id="L136" rel="#L136">136</span>
<span id="L137" rel="#L137">137</span>
<span id="L138" rel="#L138">138</span>
<span id="L139" rel="#L139">139</span>
<span id="L140" rel="#L140">140</span>
<span id="L141" rel="#L141">141</span>
<span id="L142" rel="#L142">142</span>
<span id="L143" rel="#L143">143</span>
<span id="L144" rel="#L144">144</span>
<span id="L145" rel="#L145">145</span>
<span id="L146" rel="#L146">146</span>
<span id="L147" rel="#L147">147</span>
<span id="L148" rel="#L148">148</span>
<span id="L149" rel="#L149">149</span>
<span id="L150" rel="#L150">150</span>
<span id="L151" rel="#L151">151</span>
<span id="L152" rel="#L152">152</span>
<span id="L153" rel="#L153">153</span>
<span id="L154" rel="#L154">154</span>
<span id="L155" rel="#L155">155</span>
<span id="L156" rel="#L156">156</span>
<span id="L157" rel="#L157">157</span>
<span id="L158" rel="#L158">158</span>
<span id="L159" rel="#L159">159</span>
<span id="L160" rel="#L160">160</span>
<span id="L161" rel="#L161">161</span>
<span id="L162" rel="#L162">162</span>
<span id="L163" rel="#L163">163</span>
<span id="L164" rel="#L164">164</span>
<span id="L165" rel="#L165">165</span>
<span id="L166" rel="#L166">166</span>
<span id="L167" rel="#L167">167</span>
<span id="L168" rel="#L168">168</span>
<span id="L169" rel="#L169">169</span>
<span id="L170" rel="#L170">170</span>
<span id="L171" rel="#L171">171</span>
<span id="L172" rel="#L172">172</span>
<span id="L173" rel="#L173">173</span>
<span id="L174" rel="#L174">174</span>
<span id="L175" rel="#L175">175</span>
<span id="L176" rel="#L176">176</span>
<span id="L177" rel="#L177">177</span>
<span id="L178" rel="#L178">178</span>
<span id="L179" rel="#L179">179</span>
<span id="L180" rel="#L180">180</span>
<span id="L181" rel="#L181">181</span>
<span id="L182" rel="#L182">182</span>
<span id="L183" rel="#L183">183</span>
<span id="L184" rel="#L184">184</span>
<span id="L185" rel="#L185">185</span>
<span id="L186" rel="#L186">186</span>
<span id="L187" rel="#L187">187</span>
<span id="L188" rel="#L188">188</span>
<span id="L189" rel="#L189">189</span>
<span id="L190" rel="#L190">190</span>
<span id="L191" rel="#L191">191</span>
<span id="L192" rel="#L192">192</span>
<span id="L193" rel="#L193">193</span>
<span id="L194" rel="#L194">194</span>
<span id="L195" rel="#L195">195</span>
<span id="L196" rel="#L196">196</span>
<span id="L197" rel="#L197">197</span>
<span id="L198" rel="#L198">198</span>
<span id="L199" rel="#L199">199</span>
<span id="L200" rel="#L200">200</span>
<span id="L201" rel="#L201">201</span>
<span id="L202" rel="#L202">202</span>
<span id="L203" rel="#L203">203</span>
<span id="L204" rel="#L204">204</span>
<span id="L205" rel="#L205">205</span>
<span id="L206" rel="#L206">206</span>
<span id="L207" rel="#L207">207</span>
<span id="L208" rel="#L208">208</span>
<span id="L209" rel="#L209">209</span>
<span id="L210" rel="#L210">210</span>
<span id="L211" rel="#L211">211</span>
<span id="L212" rel="#L212">212</span>
<span id="L213" rel="#L213">213</span>
<span id="L214" rel="#L214">214</span>
<span id="L215" rel="#L215">215</span>
<span id="L216" rel="#L216">216</span>
<span id="L217" rel="#L217">217</span>
<span id="L218" rel="#L218">218</span>
<span id="L219" rel="#L219">219</span>
<span id="L220" rel="#L220">220</span>
<span id="L221" rel="#L221">221</span>
<span id="L222" rel="#L222">222</span>
<span id="L223" rel="#L223">223</span>
<span id="L224" rel="#L224">224</span>
<span id="L225" rel="#L225">225</span>
<span id="L226" rel="#L226">226</span>
<span id="L227" rel="#L227">227</span>
<span id="L228" rel="#L228">228</span>
<span id="L229" rel="#L229">229</span>
<span id="L230" rel="#L230">230</span>
<span id="L231" rel="#L231">231</span>
<span id="L232" rel="#L232">232</span>
<span id="L233" rel="#L233">233</span>
<span id="L234" rel="#L234">234</span>
<span id="L235" rel="#L235">235</span>
<span id="L236" rel="#L236">236</span>
<span id="L237" rel="#L237">237</span>
<span id="L238" rel="#L238">238</span>
<span id="L239" rel="#L239">239</span>
<span id="L240" rel="#L240">240</span>
<span id="L241" rel="#L241">241</span>
<span id="L242" rel="#L242">242</span>
<span id="L243" rel="#L243">243</span>
<span id="L244" rel="#L244">244</span>
<span id="L245" rel="#L245">245</span>
<span id="L246" rel="#L246">246</span>
<span id="L247" rel="#L247">247</span>
<span id="L248" rel="#L248">248</span>
<span id="L249" rel="#L249">249</span>
<span id="L250" rel="#L250">250</span>
<span id="L251" rel="#L251">251</span>
<span id="L252" rel="#L252">252</span>
<span id="L253" rel="#L253">253</span>
<span id="L254" rel="#L254">254</span>
<span id="L255" rel="#L255">255</span>
<span id="L256" rel="#L256">256</span>
<span id="L257" rel="#L257">257</span>
<span id="L258" rel="#L258">258</span>
<span id="L259" rel="#L259">259</span>
<span id="L260" rel="#L260">260</span>
<span id="L261" rel="#L261">261</span>
<span id="L262" rel="#L262">262</span>
<span id="L263" rel="#L263">263</span>
<span id="L264" rel="#L264">264</span>
<span id="L265" rel="#L265">265</span>
<span id="L266" rel="#L266">266</span>
<span id="L267" rel="#L267">267</span>
<span id="L268" rel="#L268">268</span>
<span id="L269" rel="#L269">269</span>
<span id="L270" rel="#L270">270</span>
<span id="L271" rel="#L271">271</span>
<span id="L272" rel="#L272">272</span>
<span id="L273" rel="#L273">273</span>
<span id="L274" rel="#L274">274</span>
<span id="L275" rel="#L275">275</span>
<span id="L276" rel="#L276">276</span>
<span id="L277" rel="#L277">277</span>
<span id="L278" rel="#L278">278</span>
<span id="L279" rel="#L279">279</span>
<span id="L280" rel="#L280">280</span>
<span id="L281" rel="#L281">281</span>
<span id="L282" rel="#L282">282</span>
<span id="L283" rel="#L283">283</span>
<span id="L284" rel="#L284">284</span>
<span id="L285" rel="#L285">285</span>
<span id="L286" rel="#L286">286</span>
<span id="L287" rel="#L287">287</span>
<span id="L288" rel="#L288">288</span>
<span id="L289" rel="#L289">289</span>
<span id="L290" rel="#L290">290</span>
<span id="L291" rel="#L291">291</span>
<span id="L292" rel="#L292">292</span>
<span id="L293" rel="#L293">293</span>
<span id="L294" rel="#L294">294</span>
<span id="L295" rel="#L295">295</span>
<span id="L296" rel="#L296">296</span>
<span id="L297" rel="#L297">297</span>
<span id="L298" rel="#L298">298</span>
<span id="L299" rel="#L299">299</span>
<span id="L300" rel="#L300">300</span>
<span id="L301" rel="#L301">301</span>
<span id="L302" rel="#L302">302</span>
<span id="L303" rel="#L303">303</span>
<span id="L304" rel="#L304">304</span>
<span id="L305" rel="#L305">305</span>
<span id="L306" rel="#L306">306</span>
<span id="L307" rel="#L307">307</span>
<span id="L308" rel="#L308">308</span>
<span id="L309" rel="#L309">309</span>
<span id="L310" rel="#L310">310</span>
<span id="L311" rel="#L311">311</span>
<span id="L312" rel="#L312">312</span>
<span id="L313" rel="#L313">313</span>
<span id="L314" rel="#L314">314</span>
<span id="L315" rel="#L315">315</span>
<span id="L316" rel="#L316">316</span>
<span id="L317" rel="#L317">317</span>
<span id="L318" rel="#L318">318</span>
<span id="L319" rel="#L319">319</span>
<span id="L320" rel="#L320">320</span>
<span id="L321" rel="#L321">321</span>
<span id="L322" rel="#L322">322</span>
<span id="L323" rel="#L323">323</span>
<span id="L324" rel="#L324">324</span>
<span id="L325" rel="#L325">325</span>
<span id="L326" rel="#L326">326</span>
<span id="L327" rel="#L327">327</span>
<span id="L328" rel="#L328">328</span>
<span id="L329" rel="#L329">329</span>
<span id="L330" rel="#L330">330</span>
<span id="L331" rel="#L331">331</span>
<span id="L332" rel="#L332">332</span>
<span id="L333" rel="#L333">333</span>
<span id="L334" rel="#L334">334</span>
<span id="L335" rel="#L335">335</span>
<span id="L336" rel="#L336">336</span>
<span id="L337" rel="#L337">337</span>
<span id="L338" rel="#L338">338</span>
<span id="L339" rel="#L339">339</span>
<span id="L340" rel="#L340">340</span>
<span id="L341" rel="#L341">341</span>
<span id="L342" rel="#L342">342</span>
<span id="L343" rel="#L343">343</span>
<span id="L344" rel="#L344">344</span>
<span id="L345" rel="#L345">345</span>
<span id="L346" rel="#L346">346</span>
<span id="L347" rel="#L347">347</span>
<span id="L348" rel="#L348">348</span>
<span id="L349" rel="#L349">349</span>
<span id="L350" rel="#L350">350</span>
<span id="L351" rel="#L351">351</span>
<span id="L352" rel="#L352">352</span>
<span id="L353" rel="#L353">353</span>
<span id="L354" rel="#L354">354</span>
<span id="L355" rel="#L355">355</span>
<span id="L356" rel="#L356">356</span>
<span id="L357" rel="#L357">357</span>
<span id="L358" rel="#L358">358</span>
<span id="L359" rel="#L359">359</span>
<span id="L360" rel="#L360">360</span>
<span id="L361" rel="#L361">361</span>
<span id="L362" rel="#L362">362</span>
<span id="L363" rel="#L363">363</span>
<span id="L364" rel="#L364">364</span>
<span id="L365" rel="#L365">365</span>
<span id="L366" rel="#L366">366</span>
<span id="L367" rel="#L367">367</span>
<span id="L368" rel="#L368">368</span>
<span id="L369" rel="#L369">369</span>
<span id="L370" rel="#L370">370</span>
<span id="L371" rel="#L371">371</span>
<span id="L372" rel="#L372">372</span>
<span id="L373" rel="#L373">373</span>
<span id="L374" rel="#L374">374</span>
<span id="L375" rel="#L375">375</span>
<span id="L376" rel="#L376">376</span>
<span id="L377" rel="#L377">377</span>
<span id="L378" rel="#L378">378</span>
<span id="L379" rel="#L379">379</span>
<span id="L380" rel="#L380">380</span>
<span id="L381" rel="#L381">381</span>
<span id="L382" rel="#L382">382</span>
<span id="L383" rel="#L383">383</span>
<span id="L384" rel="#L384">384</span>
<span id="L385" rel="#L385">385</span>
<span id="L386" rel="#L386">386</span>
<span id="L387" rel="#L387">387</span>
<span id="L388" rel="#L388">388</span>
<span id="L389" rel="#L389">389</span>
<span id="L390" rel="#L390">390</span>
<span id="L391" rel="#L391">391</span>
<span id="L392" rel="#L392">392</span>
<span id="L393" rel="#L393">393</span>
<span id="L394" rel="#L394">394</span>
<span id="L395" rel="#L395">395</span>
<span id="L396" rel="#L396">396</span>
<span id="L397" rel="#L397">397</span>
<span id="L398" rel="#L398">398</span>
<span id="L399" rel="#L399">399</span>
<span id="L400" rel="#L400">400</span>
<span id="L401" rel="#L401">401</span>
<span id="L402" rel="#L402">402</span>
<span id="L403" rel="#L403">403</span>
<span id="L404" rel="#L404">404</span>
<span id="L405" rel="#L405">405</span>
<span id="L406" rel="#L406">406</span>
<span id="L407" rel="#L407">407</span>
<span id="L408" rel="#L408">408</span>
<span id="L409" rel="#L409">409</span>
<span id="L410" rel="#L410">410</span>
<span id="L411" rel="#L411">411</span>
<span id="L412" rel="#L412">412</span>
<span id="L413" rel="#L413">413</span>
<span id="L414" rel="#L414">414</span>
<span id="L415" rel="#L415">415</span>
<span id="L416" rel="#L416">416</span>
<span id="L417" rel="#L417">417</span>
<span id="L418" rel="#L418">418</span>
<span id="L419" rel="#L419">419</span>
<span id="L420" rel="#L420">420</span>
<span id="L421" rel="#L421">421</span>
<span id="L422" rel="#L422">422</span>
<span id="L423" rel="#L423">423</span>
<span id="L424" rel="#L424">424</span>
<span id="L425" rel="#L425">425</span>
<span id="L426" rel="#L426">426</span>
<span id="L427" rel="#L427">427</span>
<span id="L428" rel="#L428">428</span>
<span id="L429" rel="#L429">429</span>
<span id="L430" rel="#L430">430</span>
<span id="L431" rel="#L431">431</span>
<span id="L432" rel="#L432">432</span>
<span id="L433" rel="#L433">433</span>
<span id="L434" rel="#L434">434</span>
<span id="L435" rel="#L435">435</span>
<span id="L436" rel="#L436">436</span>
<span id="L437" rel="#L437">437</span>
<span id="L438" rel="#L438">438</span>
<span id="L439" rel="#L439">439</span>
<span id="L440" rel="#L440">440</span>
<span id="L441" rel="#L441">441</span>
<span id="L442" rel="#L442">442</span>
<span id="L443" rel="#L443">443</span>
<span id="L444" rel="#L444">444</span>
<span id="L445" rel="#L445">445</span>
<span id="L446" rel="#L446">446</span>
<span id="L447" rel="#L447">447</span>
<span id="L448" rel="#L448">448</span>
<span id="L449" rel="#L449">449</span>
<span id="L450" rel="#L450">450</span>
<span id="L451" rel="#L451">451</span>
<span id="L452" rel="#L452">452</span>
<span id="L453" rel="#L453">453</span>
<span id="L454" rel="#L454">454</span>
<span id="L455" rel="#L455">455</span>
<span id="L456" rel="#L456">456</span>
<span id="L457" rel="#L457">457</span>
<span id="L458" rel="#L458">458</span>
<span id="L459" rel="#L459">459</span>
<span id="L460" rel="#L460">460</span>
<span id="L461" rel="#L461">461</span>
<span id="L462" rel="#L462">462</span>
<span id="L463" rel="#L463">463</span>
<span id="L464" rel="#L464">464</span>
<span id="L465" rel="#L465">465</span>
<span id="L466" rel="#L466">466</span>
<span id="L467" rel="#L467">467</span>
<span id="L468" rel="#L468">468</span>
<span id="L469" rel="#L469">469</span>
<span id="L470" rel="#L470">470</span>
<span id="L471" rel="#L471">471</span>
<span id="L472" rel="#L472">472</span>
<span id="L473" rel="#L473">473</span>
<span id="L474" rel="#L474">474</span>
<span id="L475" rel="#L475">475</span>
<span id="L476" rel="#L476">476</span>
<span id="L477" rel="#L477">477</span>
<span id="L478" rel="#L478">478</span>
<span id="L479" rel="#L479">479</span>
<span id="L480" rel="#L480">480</span>
<span id="L481" rel="#L481">481</span>
<span id="L482" rel="#L482">482</span>
<span id="L483" rel="#L483">483</span>
<span id="L484" rel="#L484">484</span>
<span id="L485" rel="#L485">485</span>
<span id="L486" rel="#L486">486</span>
<span id="L487" rel="#L487">487</span>
<span id="L488" rel="#L488">488</span>
<span id="L489" rel="#L489">489</span>
<span id="L490" rel="#L490">490</span>
<span id="L491" rel="#L491">491</span>
<span id="L492" rel="#L492">492</span>
<span id="L493" rel="#L493">493</span>
<span id="L494" rel="#L494">494</span>
<span id="L495" rel="#L495">495</span>
<span id="L496" rel="#L496">496</span>
<span id="L497" rel="#L497">497</span>
<span id="L498" rel="#L498">498</span>
<span id="L499" rel="#L499">499</span>
<span id="L500" rel="#L500">500</span>
<span id="L501" rel="#L501">501</span>
<span id="L502" rel="#L502">502</span>
<span id="L503" rel="#L503">503</span>
<span id="L504" rel="#L504">504</span>
<span id="L505" rel="#L505">505</span>
<span id="L506" rel="#L506">506</span>
<span id="L507" rel="#L507">507</span>
<span id="L508" rel="#L508">508</span>
<span id="L509" rel="#L509">509</span>
<span id="L510" rel="#L510">510</span>
<span id="L511" rel="#L511">511</span>
<span id="L512" rel="#L512">512</span>
<span id="L513" rel="#L513">513</span>
<span id="L514" rel="#L514">514</span>
<span id="L515" rel="#L515">515</span>
<span id="L516" rel="#L516">516</span>
<span id="L517" rel="#L517">517</span>
<span id="L518" rel="#L518">518</span>
<span id="L519" rel="#L519">519</span>
<span id="L520" rel="#L520">520</span>
<span id="L521" rel="#L521">521</span>
<span id="L522" rel="#L522">522</span>
<span id="L523" rel="#L523">523</span>
<span id="L524" rel="#L524">524</span>
<span id="L525" rel="#L525">525</span>
<span id="L526" rel="#L526">526</span>
<span id="L527" rel="#L527">527</span>
<span id="L528" rel="#L528">528</span>
<span id="L529" rel="#L529">529</span>
<span id="L530" rel="#L530">530</span>
<span id="L531" rel="#L531">531</span>
<span id="L532" rel="#L532">532</span>
<span id="L533" rel="#L533">533</span>
<span id="L534" rel="#L534">534</span>
<span id="L535" rel="#L535">535</span>
<span id="L536" rel="#L536">536</span>
<span id="L537" rel="#L537">537</span>
<span id="L538" rel="#L538">538</span>
<span id="L539" rel="#L539">539</span>
<span id="L540" rel="#L540">540</span>
<span id="L541" rel="#L541">541</span>
<span id="L542" rel="#L542">542</span>
<span id="L543" rel="#L543">543</span>
<span id="L544" rel="#L544">544</span>
<span id="L545" rel="#L545">545</span>
<span id="L546" rel="#L546">546</span>
<span id="L547" rel="#L547">547</span>
<span id="L548" rel="#L548">548</span>
<span id="L549" rel="#L549">549</span>
<span id="L550" rel="#L550">550</span>
<span id="L551" rel="#L551">551</span>
<span id="L552" rel="#L552">552</span>
<span id="L553" rel="#L553">553</span>
<span id="L554" rel="#L554">554</span>
<span id="L555" rel="#L555">555</span>
<span id="L556" rel="#L556">556</span>
<span id="L557" rel="#L557">557</span>
<span id="L558" rel="#L558">558</span>
<span id="L559" rel="#L559">559</span>
<span id="L560" rel="#L560">560</span>
<span id="L561" rel="#L561">561</span>
<span id="L562" rel="#L562">562</span>
<span id="L563" rel="#L563">563</span>
<span id="L564" rel="#L564">564</span>
<span id="L565" rel="#L565">565</span>
<span id="L566" rel="#L566">566</span>
<span id="L567" rel="#L567">567</span>
<span id="L568" rel="#L568">568</span>
<span id="L569" rel="#L569">569</span>
<span id="L570" rel="#L570">570</span>
<span id="L571" rel="#L571">571</span>
<span id="L572" rel="#L572">572</span>
<span id="L573" rel="#L573">573</span>
<span id="L574" rel="#L574">574</span>
<span id="L575" rel="#L575">575</span>
<span id="L576" rel="#L576">576</span>
<span id="L577" rel="#L577">577</span>
<span id="L578" rel="#L578">578</span>
<span id="L579" rel="#L579">579</span>
<span id="L580" rel="#L580">580</span>
<span id="L581" rel="#L581">581</span>
<span id="L582" rel="#L582">582</span>
<span id="L583" rel="#L583">583</span>
<span id="L584" rel="#L584">584</span>
<span id="L585" rel="#L585">585</span>
<span id="L586" rel="#L586">586</span>
<span id="L587" rel="#L587">587</span>
<span id="L588" rel="#L588">588</span>
<span id="L589" rel="#L589">589</span>
<span id="L590" rel="#L590">590</span>
<span id="L591" rel="#L591">591</span>
<span id="L592" rel="#L592">592</span>
<span id="L593" rel="#L593">593</span>
<span id="L594" rel="#L594">594</span>
<span id="L595" rel="#L595">595</span>
<span id="L596" rel="#L596">596</span>
<span id="L597" rel="#L597">597</span>
<span id="L598" rel="#L598">598</span>
<span id="L599" rel="#L599">599</span>
<span id="L600" rel="#L600">600</span>
<span id="L601" rel="#L601">601</span>
<span id="L602" rel="#L602">602</span>
<span id="L603" rel="#L603">603</span>
<span id="L604" rel="#L604">604</span>
<span id="L605" rel="#L605">605</span>
<span id="L606" rel="#L606">606</span>
<span id="L607" rel="#L607">607</span>
<span id="L608" rel="#L608">608</span>
<span id="L609" rel="#L609">609</span>
<span id="L610" rel="#L610">610</span>
<span id="L611" rel="#L611">611</span>
<span id="L612" rel="#L612">612</span>
<span id="L613" rel="#L613">613</span>
<span id="L614" rel="#L614">614</span>
<span id="L615" rel="#L615">615</span>
<span id="L616" rel="#L616">616</span>
<span id="L617" rel="#L617">617</span>
<span id="L618" rel="#L618">618</span>
<span id="L619" rel="#L619">619</span>
<span id="L620" rel="#L620">620</span>
<span id="L621" rel="#L621">621</span>
<span id="L622" rel="#L622">622</span>
<span id="L623" rel="#L623">623</span>
<span id="L624" rel="#L624">624</span>
<span id="L625" rel="#L625">625</span>
<span id="L626" rel="#L626">626</span>
<span id="L627" rel="#L627">627</span>
<span id="L628" rel="#L628">628</span>
<span id="L629" rel="#L629">629</span>
<span id="L630" rel="#L630">630</span>
<span id="L631" rel="#L631">631</span>
<span id="L632" rel="#L632">632</span>

           </td>
           <td class="blob-line-code"><div class="code-body highlight"><pre><div class='line' id='LC1'><span class="cp">#include &quot;Rivet/Analysis.hh&quot;</span></div><div class='line' id='LC2'><span class="cp">#include &quot;Rivet/RivetAIDA.hh&quot;</span></div><div class='line' id='LC3'><span class="cp">#include &quot;Rivet/Tools/Logging.hh&quot;</span></div><div class='line' id='LC4'><span class="cp">#include &quot;Rivet/Projections/FinalState.hh&quot;</span></div><div class='line' id='LC5'><span class="cp">#include &quot;Rivet/Projections/FastJets.hh&quot;</span></div><div class='line' id='LC6'><span class="cp">#include &quot;Rivet/Projections/ChargedFinalState.hh&quot;</span></div><div class='line' id='LC7'><span class="cp">#include &quot;Rivet/Projections/InvMassFinalState.hh&quot;</span></div><div class='line' id='LC8'><span class="cp">#include &quot;Rivet/Tools/ParticleIdUtils.hh&quot;</span></div><div class='line' id='LC9'><span class="cp">#include &quot;Rivet/Math/Vector4.hh&quot;</span></div><div class='line' id='LC10'><br/></div><div class='line' id='LC11'><span class="c1">//#define DebugLog </span></div><div class='line' id='LC12'><br/></div><div class='line' id='LC13'><span class="k">namespace</span> <span class="n">Rivet</span></div><div class='line' id='LC14'><span class="p">{</span></div><div class='line' id='LC15'><br/></div><div class='line' id='LC16'>&nbsp;&nbsp;<span class="k">class</span> <span class="nc">CMS_SMP_12_017</span><span class="o">:</span><span class="k">public</span> <span class="n">Analysis</span></div><div class='line' id='LC17'>&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC18'>&nbsp;&nbsp;<span class="k">public</span><span class="o">:</span></div><div class='line' id='LC19'><br/></div><div class='line' id='LC20'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">/// Constructor</span></div><div class='line' id='LC21'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">CMS_SMP_12_017</span><span class="p">()</span><span class="o">:</span><span class="n">Analysis</span><span class="p">(</span><span class="s">&quot;CMS_SMP_12_017&quot;</span><span class="p">)</span></div><div class='line' id='LC22'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC23'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">setNeedsCrossSection</span><span class="p">(</span><span class="nb">true</span><span class="p">);</span></div><div class='line' id='LC24'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC25'><br/></div><div class='line' id='LC26'>&nbsp;&nbsp;<span class="k">public</span><span class="o">:</span></div><div class='line' id='LC27'><br/></div><div class='line' id='LC28'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">/// Book histograms and initialise projections before the run</span></div><div class='line' id='LC29'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">void</span> <span class="n">init</span><span class="p">()</span></div><div class='line' id='LC30'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC31'><br/></div><div class='line' id='LC32'><span class="cp">#ifdef DebugLog</span></div><div class='line' id='LC33'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// set optionally the verbosity for the internal Rivet message system</span></div><div class='line' id='LC34'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">getLog</span><span class="p">().</span><span class="n">setLevel</span><span class="p">(</span><span class="mi">0</span><span class="p">);</span></div><div class='line' id='LC35'><span class="cp">#endif</span></div><div class='line' id='LC36'><br/></div><div class='line' id='LC37'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">const</span> <span class="n">FinalState</span> <span class="n">fs</span><span class="p">;</span></div><div class='line' id='LC38'><br/></div><div class='line' id='LC39'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">addProjection</span><span class="p">(</span><span class="n">fs</span><span class="p">,</span> <span class="s">&quot;FS&quot;</span><span class="p">);</span></div><div class='line' id='LC40'><br/></div><div class='line' id='LC41'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">SumW2</span> <span class="o">=</span> <span class="mf">0.</span><span class="p">;</span></div><div class='line' id='LC42'><br/></div><div class='line' id='LC43'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_excmult_jets_ele</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC44'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_incmult_jets_ele</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">2</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC45'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_leading_jet_pt_ele</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">3</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC46'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_second_jet_pt_ele</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">4</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC47'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_third_jet_pt_ele</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">5</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC48'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_fourth_jet_pt_ele</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">6</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC49'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_leading_jet_eta_ele</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">7</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC50'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_second_jet_eta_ele</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">8</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC51'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_third_jet_eta_ele</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">9</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC52'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_fourth_jet_eta_ele</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">10</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC53'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_ht1_ele</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">11</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC54'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_ht2_ele</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">12</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC55'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_ht3_ele</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">13</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC56'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_ht4_ele</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">14</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC57'><br/></div><div class='line' id='LC58'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_excmult_jets_muo</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">15</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC59'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_incmult_jets_muo</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">16</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC60'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_leading_jet_pt_muo</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">17</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC61'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_second_jet_pt_muo</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">18</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC62'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_third_jet_pt_muo</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">19</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC63'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_fourth_jet_pt_muo</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">20</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC64'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_leading_jet_eta_muo</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">21</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC65'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_second_jet_eta_muo</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">22</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC66'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_third_jet_eta_muo</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">23</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC67'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_fourth_jet_eta_muo</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">24</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC68'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_ht1_muo</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">25</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC69'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_ht2_muo</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">26</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC70'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_ht3_muo</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">27</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC71'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_ht4_muo</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">28</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC72'><br/></div><div class='line' id='LC73'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_excmult_jets_tot</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">29</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC74'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_incmult_jets_tot</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">30</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC75'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_leading_jet_pt_tot</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">31</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC76'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_second_jet_pt_tot</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">32</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC77'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_third_jet_pt_tot</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">33</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC78'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_fourth_jet_pt_tot</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">34</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC79'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_leading_jet_eta_tot</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">35</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC80'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_second_jet_eta_tot</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">36</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC81'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_third_jet_eta_tot</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">37</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC82'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_fourth_jet_eta_tot</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">38</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC83'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_ht1_tot</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">39</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC84'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_ht2_tot</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">40</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC85'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_ht3_tot</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">41</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC86'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_ht4_tot</span> <span class="o">=</span> <span class="n">bookHistogram1D</span><span class="p">(</span><span class="mi">42</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></div><div class='line' id='LC87'><br/></div><div class='line' id='LC88'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC89'><br/></div><div class='line' id='LC90'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">/// Perform the per-event analysis</span></div><div class='line' id='LC91'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">void</span> <span class="n">analyze</span><span class="p">(</span><span class="k">const</span> <span class="n">Event</span> <span class="o">&amp;</span> <span class="n">event</span><span class="p">)</span></div><div class='line' id='LC92'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC93'><br/></div><div class='line' id='LC94'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">const</span> <span class="kt">double</span> <span class="n">weight</span> <span class="o">=</span> <span class="n">event</span><span class="p">.</span><span class="n">weight</span><span class="p">();</span></div><div class='line' id='LC95'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">std</span><span class="o">::</span><span class="n">vector</span><span class="o">&lt;</span><span class="n">fastjet</span><span class="o">::</span><span class="n">PseudoJet</span><span class="o">&gt;</span> <span class="n">vecs</span><span class="p">;</span></div><div class='line' id='LC96'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">bool</span> <span class="n">Z_ele</span> <span class="o">=</span> <span class="nb">false</span><span class="p">;</span></div><div class='line' id='LC97'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">bool</span> <span class="n">Z_muon</span> <span class="o">=</span> <span class="nb">false</span><span class="p">;</span></div><div class='line' id='LC98'><br/></div><div class='line' id='LC99'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">SumW2</span> <span class="o">+=</span> <span class="n">weight</span> <span class="o">*</span> <span class="n">weight</span><span class="p">;</span></div><div class='line' id='LC100'><br/></div><div class='line' id='LC101'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">vector</span><span class="o">&lt;</span><span class="kt">unsigned</span> <span class="kt">int</span><span class="o">&gt;</span> <span class="n">ele_photons</span><span class="p">;</span></div><div class='line' id='LC102'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">vector</span><span class="o">&lt;</span><span class="kt">unsigned</span> <span class="kt">int</span><span class="o">&gt;</span> <span class="n">pos_photons</span><span class="p">;</span></div><div class='line' id='LC103'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">vector</span><span class="o">&lt;</span><span class="kt">unsigned</span> <span class="kt">int</span><span class="o">&gt;</span> <span class="n">muon_photons</span><span class="p">;</span></div><div class='line' id='LC104'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">vector</span><span class="o">&lt;</span><span class="kt">unsigned</span> <span class="kt">int</span><span class="o">&gt;</span> <span class="n">antimuon_photons</span><span class="p">;</span></div><div class='line' id='LC105'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">vector</span><span class="o">&lt;</span><span class="kt">unsigned</span> <span class="kt">int</span><span class="o">&gt;</span> <span class="n">photons</span><span class="p">;</span></div><div class='line' id='LC106'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">vector</span><span class="o">&lt;</span><span class="kt">unsigned</span> <span class="kt">int</span><span class="o">&gt;</span> <span class="n">neutrinos</span><span class="p">;</span></div><div class='line' id='LC107'><br/></div><div class='line' id='LC108'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">struct</span> <span class="n">pt_and_particles</span></div><div class='line' id='LC109'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC110'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">FourMomentum</span> <span class="n">p_part</span><span class="p">;</span></div><div class='line' id='LC111'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">vector</span><span class="o">&lt;</span><span class="kt">unsigned</span> <span class="kt">int</span><span class="o">&gt;</span><span class="n">lepton_photon</span><span class="p">;</span></div><div class='line' id='LC112'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">};</span></div><div class='line' id='LC113'><br/></div><div class='line' id='LC114'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// structure containing the lepton_photon vector index of the lepton candidate daughter of Z and dressed lepton kinematics</span></div><div class='line' id='LC115'><br/></div><div class='line' id='LC116'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">struct</span> <span class="n">pt_and_particles</span> <span class="n">ele_dres</span><span class="p">;</span></div><div class='line' id='LC117'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">struct</span> <span class="n">pt_and_particles</span> <span class="n">pos_dres</span><span class="p">;</span></div><div class='line' id='LC118'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">struct</span> <span class="n">pt_and_particles</span> <span class="n">muon_dres</span><span class="p">;</span></div><div class='line' id='LC119'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">struct</span> <span class="n">pt_and_particles</span> <span class="n">antimuon_dres</span><span class="p">;</span></div><div class='line' id='LC120'><br/></div><div class='line' id='LC121'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">vector</span><span class="o">&lt;</span><span class="n">Particle</span><span class="o">&gt;</span> <span class="n">part_jets</span><span class="p">;</span></div><div class='line' id='LC122'><br/></div><div class='line' id='LC123'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">const</span> <span class="n">FinalState</span> <span class="o">&amp;</span> <span class="n">fs</span> <span class="o">=</span> <span class="n">applyProjection</span><span class="o">&lt;</span> <span class="n">FinalState</span><span class="o">&gt;</span><span class="p">(</span><span class="n">event</span><span class="p">,</span> <span class="s">&quot;FS&quot;</span><span class="p">);</span></div><div class='line' id='LC124'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">fs</span><span class="p">.</span><span class="n">particlesByPt</span><span class="p">();</span></div><div class='line' id='LC125'><br/></div><div class='line' id='LC126'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// copy of the list of final state partcles</span></div><div class='line' id='LC127'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">part_jets</span> <span class="o">=</span> <span class="n">fs</span><span class="p">.</span><span class="n">particles</span><span class="p">();</span></div><div class='line' id='LC128'><br/></div><div class='line' id='LC129'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">for</span> <span class="p">(</span><span class="kt">unsigned</span> <span class="kt">int</span> <span class="n">i</span> <span class="o">=</span> <span class="mi">0</span><span class="p">;</span> <span class="n">i</span> <span class="o">&lt;</span> <span class="n">fs</span><span class="p">.</span><span class="n">particles</span><span class="p">().</span><span class="n">size</span><span class="p">();</span> <span class="n">i</span><span class="o">++</span><span class="p">)</span></div><div class='line' id='LC130'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC131'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">fs</span><span class="p">.</span><span class="n">particles</span><span class="p">().</span><span class="n">at</span><span class="p">(</span><span class="n">i</span><span class="p">).</span><span class="n">pdgId</span><span class="p">()</span> <span class="o">==</span> <span class="n">PHOTON</span><span class="p">)</span></div><div class='line' id='LC132'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">photons</span><span class="p">.</span><span class="n">push_back</span><span class="p">(</span><span class="n">i</span><span class="p">);</span></div><div class='line' id='LC133'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">fabs</span><span class="p">(</span><span class="n">fs</span><span class="p">.</span><span class="n">particles</span><span class="p">().</span><span class="n">at</span><span class="p">(</span><span class="n">i</span><span class="p">).</span><span class="n">pdgId</span><span class="p">())</span> <span class="o">==</span> <span class="mi">12</span> <span class="o">||</span></div><div class='line' id='LC134'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">fabs</span><span class="p">(</span><span class="n">fs</span><span class="p">.</span><span class="n">particles</span><span class="p">().</span><span class="n">at</span><span class="p">(</span><span class="n">i</span><span class="p">).</span><span class="n">pdgId</span><span class="p">())</span> <span class="o">==</span> <span class="mi">14</span> <span class="o">||</span></div><div class='line' id='LC135'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">fabs</span><span class="p">(</span><span class="n">fs</span><span class="p">.</span><span class="n">particles</span><span class="p">().</span><span class="n">at</span><span class="p">(</span><span class="n">i</span><span class="p">).</span><span class="n">pdgId</span><span class="p">())</span> <span class="o">==</span> <span class="mi">16</span><span class="p">)</span></div><div class='line' id='LC136'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">neutrinos</span><span class="p">.</span><span class="n">push_back</span><span class="p">(</span><span class="n">i</span><span class="p">);</span></div><div class='line' id='LC137'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC138'><br/></div><div class='line' id='LC139'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">for</span> <span class="p">(</span><span class="kt">unsigned</span> <span class="kt">int</span> <span class="n">i</span> <span class="o">=</span> <span class="mi">0</span><span class="p">;</span> <span class="n">i</span> <span class="o">&lt;</span> <span class="n">fs</span><span class="p">.</span><span class="n">particles</span><span class="p">().</span><span class="n">size</span><span class="p">();</span> <span class="n">i</span><span class="o">++</span><span class="p">)</span></div><div class='line' id='LC140'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC141'><br/></div><div class='line' id='LC142'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">eta</span> <span class="o">=</span> <span class="n">fs</span><span class="p">.</span><span class="n">particles</span><span class="p">().</span><span class="n">at</span><span class="p">(</span><span class="n">i</span><span class="p">).</span><span class="n">momentum</span><span class="p">().</span><span class="n">eta</span><span class="p">();</span></div><div class='line' id='LC143'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">phi</span> <span class="o">=</span> <span class="n">fs</span><span class="p">.</span><span class="n">particles</span><span class="p">().</span><span class="n">at</span><span class="p">(</span><span class="n">i</span><span class="p">).</span><span class="n">momentum</span><span class="p">().</span><span class="n">phi</span><span class="p">();</span></div><div class='line' id='LC144'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">int</span> <span class="n">Id</span> <span class="o">=</span> <span class="n">fs</span><span class="p">.</span><span class="n">particles</span><span class="p">().</span><span class="n">at</span><span class="p">(</span><span class="n">i</span><span class="p">).</span><span class="n">pdgId</span><span class="p">();</span></div><div class='line' id='LC145'><br/></div><div class='line' id='LC146'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// electron selection</span></div><div class='line' id='LC147'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">Id</span> <span class="o">==</span> <span class="n">ELECTRON</span> <span class="o">&amp;&amp;</span> <span class="n">fabs</span><span class="p">(</span><span class="n">eta</span><span class="p">)</span> <span class="o">&lt;</span> <span class="mf">2.4</span><span class="p">)</span></div><div class='line' id='LC148'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC149'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">FourMomentum</span> <span class="n">ele_p</span><span class="p">(</span><span class="n">fs</span><span class="p">.</span><span class="n">particles</span><span class="p">().</span><span class="n">at</span><span class="p">(</span><span class="n">i</span><span class="p">).</span><span class="n">momentum</span><span class="p">());</span></div><div class='line' id='LC150'><br/></div><div class='line' id='LC151'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">ele_photons</span><span class="p">.</span><span class="n">push_back</span><span class="p">(</span><span class="n">i</span><span class="p">);</span></div><div class='line' id='LC152'><br/></div><div class='line' id='LC153'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">//lepton dressing</span></div><div class='line' id='LC154'><br/></div><div class='line' id='LC155'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">for</span> <span class="p">(</span><span class="kt">unsigned</span> <span class="kt">int</span> <span class="n">j</span> <span class="o">=</span> <span class="mi">0</span><span class="p">;</span> <span class="n">j</span> <span class="o">&lt;</span> <span class="n">photons</span><span class="p">.</span><span class="n">size</span><span class="p">();</span> <span class="n">j</span><span class="o">++</span><span class="p">)</span></div><div class='line' id='LC156'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC157'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">eta_photon</span> <span class="o">=</span> <span class="n">fs</span><span class="p">.</span><span class="n">particles</span><span class="p">().</span><span class="n">at</span><span class="p">(</span><span class="n">photons</span><span class="p">.</span><span class="n">at</span><span class="p">(</span><span class="n">j</span><span class="p">)).</span><span class="n">momentum</span><span class="p">().</span><span class="n">eta</span><span class="p">();</span></div><div class='line' id='LC158'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">phi_photon</span> <span class="o">=</span> <span class="n">fs</span><span class="p">.</span><span class="n">particles</span><span class="p">().</span><span class="n">at</span><span class="p">(</span><span class="n">photons</span><span class="p">.</span><span class="n">at</span><span class="p">(</span><span class="n">j</span><span class="p">)).</span><span class="n">momentum</span><span class="p">().</span><span class="n">phi</span><span class="p">();</span></div><div class='line' id='LC159'><br/></div><div class='line' id='LC160'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">deltaR</span> <span class="o">=</span> <span class="n">Rivet</span><span class="o">::</span><span class="n">deltaR</span><span class="p">(</span><span class="n">eta_photon</span><span class="p">,</span><span class="n">phi_photon</span><span class="p">,</span><span class="n">eta</span><span class="p">,</span><span class="n">phi</span><span class="p">);</span></div><div class='line' id='LC161'><br/></div><div class='line' id='LC162'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">deltaR</span> <span class="o">&lt;</span> <span class="mf">0.1</span><span class="p">)</span></div><div class='line' id='LC163'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC164'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">ele_photons</span><span class="p">.</span><span class="n">push_back</span><span class="p">(</span><span class="n">photons</span><span class="p">.</span><span class="n">at</span><span class="p">(</span><span class="n">j</span><span class="p">));</span></div><div class='line' id='LC165'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">ele_p</span> <span class="o">+=</span> <span class="n">fs</span><span class="p">.</span><span class="n">particles</span><span class="p">().</span><span class="n">at</span><span class="p">(</span><span class="n">photons</span><span class="p">.</span><span class="n">at</span><span class="p">(</span><span class="n">j</span><span class="p">)).</span><span class="n">momentum</span><span class="p">();</span></div><div class='line' id='LC166'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC167'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC168'><br/></div><div class='line' id='LC169'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// keep dressed electron if above threshold</span></div><div class='line' id='LC170'><br/></div><div class='line' id='LC171'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">ele_dres</span><span class="p">.</span><span class="n">lepton_photon</span><span class="p">.</span><span class="n">empty</span><span class="p">()</span> <span class="o">&amp;&amp;</span> <span class="n">ele_p</span><span class="p">.</span><span class="n">pT</span><span class="p">()</span> <span class="o">&gt;</span> <span class="mi">20</span><span class="p">)</span></div><div class='line' id='LC172'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC173'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">ele_dres</span><span class="p">.</span><span class="n">p_part</span> <span class="o">=</span> <span class="n">ele_p</span><span class="p">;</span></div><div class='line' id='LC174'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">ele_dres</span><span class="p">.</span><span class="n">lepton_photon</span> <span class="o">=</span> <span class="n">ele_photons</span><span class="p">;</span></div><div class='line' id='LC175'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC176'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">else</span></div><div class='line' id='LC177'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC178'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">ele_p</span><span class="p">.</span><span class="n">pT</span><span class="p">()</span> <span class="o">&gt;</span> <span class="n">ele_dres</span><span class="p">.</span><span class="n">p_part</span><span class="p">.</span><span class="n">pT</span><span class="p">()</span> <span class="o">&amp;&amp;</span> <span class="n">ele_p</span><span class="p">.</span><span class="n">pT</span><span class="p">()</span> <span class="o">&gt;</span> <span class="mi">20</span><span class="p">)</span></div><div class='line' id='LC179'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC180'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">ele_dres</span><span class="p">.</span><span class="n">p_part</span> <span class="o">=</span> <span class="n">ele_p</span><span class="p">;</span></div><div class='line' id='LC181'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">ele_dres</span><span class="p">.</span><span class="n">lepton_photon</span> <span class="o">=</span> <span class="n">ele_photons</span><span class="p">;</span></div><div class='line' id='LC182'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC183'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC184'><br/></div><div class='line' id='LC185'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC186'><br/></div><div class='line' id='LC187'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// positron selection</span></div><div class='line' id='LC188'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">Id</span> <span class="o">==</span> <span class="n">POSITRON</span> <span class="o">&amp;&amp;</span> <span class="n">fabs</span><span class="p">(</span><span class="n">eta</span><span class="p">)</span> <span class="o">&lt;</span> <span class="mf">2.4</span><span class="p">)</span></div><div class='line' id='LC189'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC190'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">FourMomentum</span> <span class="n">pos_p</span><span class="p">(</span><span class="n">fs</span><span class="p">.</span><span class="n">particles</span><span class="p">().</span><span class="n">at</span><span class="p">(</span><span class="n">i</span><span class="p">).</span><span class="n">momentum</span><span class="p">());</span></div><div class='line' id='LC191'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">pos_photons</span><span class="p">.</span><span class="n">push_back</span><span class="p">(</span><span class="n">i</span><span class="p">);</span></div><div class='line' id='LC192'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">for</span> <span class="p">(</span><span class="kt">unsigned</span> <span class="kt">int</span> <span class="n">j</span> <span class="o">=</span> <span class="mi">0</span><span class="p">;</span> <span class="n">j</span> <span class="o">&lt;</span> <span class="n">photons</span><span class="p">.</span><span class="n">size</span><span class="p">();</span> <span class="n">j</span><span class="o">++</span><span class="p">)</span></div><div class='line' id='LC193'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC194'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">eta_photon</span> <span class="o">=</span> <span class="n">fs</span><span class="p">.</span><span class="n">particles</span><span class="p">().</span><span class="n">at</span><span class="p">(</span><span class="n">photons</span><span class="p">.</span><span class="n">at</span><span class="p">(</span><span class="n">j</span><span class="p">)).</span><span class="n">momentum</span><span class="p">().</span><span class="n">eta</span><span class="p">();</span></div><div class='line' id='LC195'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">phi_photon</span> <span class="o">=</span> <span class="n">fs</span><span class="p">.</span><span class="n">particles</span><span class="p">().</span><span class="n">at</span><span class="p">(</span><span class="n">photons</span><span class="p">.</span><span class="n">at</span><span class="p">(</span><span class="n">j</span><span class="p">)).</span><span class="n">momentum</span><span class="p">().</span><span class="n">phi</span><span class="p">();</span></div><div class='line' id='LC196'><br/></div><div class='line' id='LC197'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">deltaR</span> <span class="o">=</span> <span class="n">Rivet</span><span class="o">::</span><span class="n">deltaR</span><span class="p">(</span><span class="n">eta_photon</span><span class="p">,</span><span class="n">phi_photon</span><span class="p">,</span><span class="n">eta</span><span class="p">,</span><span class="n">phi</span><span class="p">);</span></div><div class='line' id='LC198'><br/></div><div class='line' id='LC199'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">deltaR</span> <span class="o">&lt;</span> <span class="mf">0.1</span><span class="p">)</span></div><div class='line' id='LC200'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC201'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">pos_photons</span><span class="p">.</span><span class="n">push_back</span><span class="p">(</span><span class="n">photons</span><span class="p">.</span><span class="n">at</span><span class="p">(</span><span class="n">j</span><span class="p">));</span></div><div class='line' id='LC202'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">pos_p</span> <span class="o">+=</span> <span class="n">fs</span><span class="p">.</span><span class="n">particles</span><span class="p">().</span><span class="n">at</span><span class="p">(</span><span class="n">photons</span><span class="p">.</span><span class="n">at</span><span class="p">(</span><span class="n">j</span><span class="p">)).</span><span class="n">momentum</span><span class="p">();</span></div><div class='line' id='LC203'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC204'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC205'><br/></div><div class='line' id='LC206'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">pos_dres</span><span class="p">.</span><span class="n">lepton_photon</span><span class="p">.</span><span class="n">empty</span><span class="p">()</span> <span class="o">&amp;&amp;</span> <span class="n">pos_p</span><span class="p">.</span><span class="n">pT</span><span class="p">()</span> <span class="o">&gt;</span> <span class="mi">20</span><span class="p">)</span></div><div class='line' id='LC207'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC208'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">pos_dres</span><span class="p">.</span><span class="n">p_part</span> <span class="o">=</span> <span class="n">pos_p</span><span class="p">;</span></div><div class='line' id='LC209'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">pos_dres</span><span class="p">.</span><span class="n">lepton_photon</span> <span class="o">=</span> <span class="n">pos_photons</span><span class="p">;</span></div><div class='line' id='LC210'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC211'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">else</span></div><div class='line' id='LC212'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC213'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">pos_p</span><span class="p">.</span><span class="n">pT</span><span class="p">()</span> <span class="o">&gt;</span> <span class="n">pos_dres</span><span class="p">.</span><span class="n">p_part</span><span class="p">.</span><span class="n">pT</span><span class="p">()</span> <span class="o">&amp;&amp;</span> <span class="n">pos_p</span><span class="p">.</span><span class="n">pT</span><span class="p">()</span> <span class="o">&gt;</span> <span class="mi">20</span><span class="p">)</span></div><div class='line' id='LC214'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC215'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">pos_dres</span><span class="p">.</span><span class="n">p_part</span> <span class="o">=</span> <span class="n">pos_p</span><span class="p">;</span></div><div class='line' id='LC216'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">pos_dres</span><span class="p">.</span><span class="n">lepton_photon</span> <span class="o">=</span> <span class="n">pos_photons</span><span class="p">;</span></div><div class='line' id='LC217'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC218'><br/></div><div class='line' id='LC219'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC220'><br/></div><div class='line' id='LC221'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC222'><br/></div><div class='line' id='LC223'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// muon selection</span></div><div class='line' id='LC224'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">Id</span> <span class="o">==</span> <span class="n">MUON</span> <span class="o">&amp;&amp;</span> <span class="n">fabs</span><span class="p">(</span><span class="n">eta</span><span class="p">)</span> <span class="o">&lt;</span> <span class="mf">2.4</span><span class="p">)</span></div><div class='line' id='LC225'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC226'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">FourMomentum</span> <span class="n">muon_p</span><span class="p">(</span><span class="n">fs</span><span class="p">.</span><span class="n">particles</span><span class="p">().</span><span class="n">at</span><span class="p">(</span><span class="n">i</span><span class="p">).</span><span class="n">momentum</span><span class="p">());</span></div><div class='line' id='LC227'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">muon_photons</span><span class="p">.</span><span class="n">push_back</span><span class="p">(</span><span class="n">i</span><span class="p">);</span></div><div class='line' id='LC228'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">for</span> <span class="p">(</span><span class="kt">unsigned</span> <span class="kt">int</span> <span class="n">j</span> <span class="o">=</span> <span class="mi">0</span><span class="p">;</span> <span class="n">j</span> <span class="o">&lt;</span> <span class="n">photons</span><span class="p">.</span><span class="n">size</span><span class="p">();</span> <span class="n">j</span><span class="o">++</span><span class="p">)</span></div><div class='line' id='LC229'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC230'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">eta_photon</span> <span class="o">=</span> <span class="n">fs</span><span class="p">.</span><span class="n">particles</span><span class="p">().</span><span class="n">at</span><span class="p">(</span><span class="n">photons</span><span class="p">.</span><span class="n">at</span><span class="p">(</span><span class="n">j</span><span class="p">)).</span><span class="n">momentum</span><span class="p">().</span><span class="n">eta</span><span class="p">();</span></div><div class='line' id='LC231'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">phi_photon</span> <span class="o">=</span> <span class="n">fs</span><span class="p">.</span><span class="n">particles</span><span class="p">().</span><span class="n">at</span><span class="p">(</span><span class="n">photons</span><span class="p">.</span><span class="n">at</span><span class="p">(</span><span class="n">j</span><span class="p">)).</span><span class="n">momentum</span><span class="p">().</span><span class="n">phi</span><span class="p">();</span></div><div class='line' id='LC232'><br/></div><div class='line' id='LC233'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">deltaR</span> <span class="o">=</span> <span class="n">Rivet</span><span class="o">::</span><span class="n">deltaR</span><span class="p">(</span><span class="n">eta_photon</span><span class="p">,</span><span class="n">phi_photon</span><span class="p">,</span><span class="n">eta</span><span class="p">,</span><span class="n">phi</span><span class="p">);</span></div><div class='line' id='LC234'><br/></div><div class='line' id='LC235'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">deltaR</span> <span class="o">&lt;</span> <span class="mf">0.1</span><span class="p">)</span></div><div class='line' id='LC236'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC237'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">muon_photons</span><span class="p">.</span><span class="n">push_back</span><span class="p">(</span><span class="n">photons</span><span class="p">.</span><span class="n">at</span><span class="p">(</span><span class="n">j</span><span class="p">));</span></div><div class='line' id='LC238'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">muon_p</span> <span class="o">+=</span> <span class="n">fs</span><span class="p">.</span><span class="n">particles</span><span class="p">().</span><span class="n">at</span><span class="p">(</span><span class="n">photons</span><span class="p">.</span><span class="n">at</span><span class="p">(</span><span class="n">j</span><span class="p">)).</span><span class="n">momentum</span><span class="p">();</span></div><div class='line' id='LC239'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC240'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC241'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC242'><br/></div><div class='line' id='LC243'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">muon_dres</span><span class="p">.</span><span class="n">lepton_photon</span><span class="p">.</span><span class="n">empty</span><span class="p">()</span> <span class="o">&amp;&amp;</span> <span class="n">muon_p</span><span class="p">.</span><span class="n">pT</span><span class="p">()</span> <span class="o">&gt;</span> <span class="mi">20</span><span class="p">)</span></div><div class='line' id='LC244'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC245'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">muon_dres</span><span class="p">.</span><span class="n">p_part</span> <span class="o">=</span> <span class="n">muon_p</span><span class="p">;</span></div><div class='line' id='LC246'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">muon_dres</span><span class="p">.</span><span class="n">lepton_photon</span> <span class="o">=</span> <span class="n">muon_photons</span><span class="p">;</span></div><div class='line' id='LC247'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC248'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">else</span></div><div class='line' id='LC249'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC250'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">muon_p</span><span class="p">.</span><span class="n">pT</span><span class="p">()</span> <span class="o">&gt;</span> <span class="n">muon_dres</span><span class="p">.</span><span class="n">p_part</span><span class="p">.</span><span class="n">pT</span><span class="p">()</span></div><div class='line' id='LC251'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="o">&amp;&amp;</span> <span class="n">muon_p</span><span class="p">.</span><span class="n">pT</span><span class="p">()</span> <span class="o">&gt;</span> <span class="mi">20</span><span class="p">)</span></div><div class='line' id='LC252'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC253'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">muon_dres</span><span class="p">.</span><span class="n">p_part</span> <span class="o">=</span> <span class="n">muon_p</span><span class="p">;</span></div><div class='line' id='LC254'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">muon_dres</span><span class="p">.</span><span class="n">lepton_photon</span> <span class="o">=</span> <span class="n">muon_photons</span><span class="p">;</span></div><div class='line' id='LC255'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC256'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC257'><br/></div><div class='line' id='LC258'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC259'><br/></div><div class='line' id='LC260'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// antimuon selection</span></div><div class='line' id='LC261'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">Id</span> <span class="o">==</span> <span class="n">ANTIMUON</span> <span class="o">&amp;&amp;</span> <span class="n">fabs</span><span class="p">(</span><span class="n">eta</span><span class="p">)</span> <span class="o">&lt;</span> <span class="mf">2.4</span><span class="p">)</span></div><div class='line' id='LC262'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC263'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">FourMomentum</span> <span class="n">muon_p</span><span class="p">(</span><span class="n">fs</span><span class="p">.</span><span class="n">particles</span><span class="p">().</span><span class="n">at</span><span class="p">(</span><span class="n">i</span><span class="p">).</span><span class="n">momentum</span><span class="p">());</span></div><div class='line' id='LC264'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">antimuon_photons</span><span class="p">.</span><span class="n">push_back</span><span class="p">(</span><span class="n">i</span><span class="p">);</span></div><div class='line' id='LC265'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">for</span> <span class="p">(</span><span class="kt">unsigned</span> <span class="kt">int</span> <span class="n">j</span> <span class="o">=</span> <span class="mi">0</span><span class="p">;</span> <span class="n">j</span> <span class="o">&lt;</span> <span class="n">photons</span><span class="p">.</span><span class="n">size</span><span class="p">();</span> <span class="n">j</span><span class="o">++</span><span class="p">)</span></div><div class='line' id='LC266'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC267'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">eta_photon</span> <span class="o">=</span> <span class="n">fs</span><span class="p">.</span><span class="n">particles</span><span class="p">().</span><span class="n">at</span><span class="p">(</span><span class="n">photons</span><span class="p">.</span><span class="n">at</span><span class="p">(</span><span class="n">j</span><span class="p">)).</span><span class="n">momentum</span><span class="p">().</span><span class="n">eta</span><span class="p">();</span></div><div class='line' id='LC268'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">phi_photon</span> <span class="o">=</span> <span class="n">fs</span><span class="p">.</span><span class="n">particles</span><span class="p">().</span><span class="n">at</span><span class="p">(</span><span class="n">photons</span><span class="p">.</span><span class="n">at</span><span class="p">(</span><span class="n">j</span><span class="p">)).</span><span class="n">momentum</span><span class="p">().</span><span class="n">phi</span><span class="p">();</span></div><div class='line' id='LC269'><br/></div><div class='line' id='LC270'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">deltaR</span> <span class="o">=</span> <span class="n">Rivet</span><span class="o">::</span><span class="n">deltaR</span><span class="p">(</span><span class="n">eta_photon</span><span class="p">,</span><span class="n">phi_photon</span><span class="p">,</span><span class="n">eta</span><span class="p">,</span><span class="n">phi</span><span class="p">);</span></div><div class='line' id='LC271'><br/></div><div class='line' id='LC272'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">deltaR</span> <span class="o">&lt;</span> <span class="mf">0.1</span><span class="p">)</span></div><div class='line' id='LC273'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC274'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">antimuon_photons</span><span class="p">.</span><span class="n">push_back</span><span class="p">(</span><span class="n">photons</span><span class="p">.</span><span class="n">at</span><span class="p">(</span><span class="n">j</span><span class="p">));</span></div><div class='line' id='LC275'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">muon_p</span> <span class="o">+=</span> <span class="n">fs</span><span class="p">.</span><span class="n">particles</span><span class="p">().</span><span class="n">at</span><span class="p">(</span><span class="n">photons</span><span class="p">.</span><span class="n">at</span><span class="p">(</span><span class="n">j</span><span class="p">)).</span><span class="n">momentum</span><span class="p">();</span></div><div class='line' id='LC276'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC277'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC278'><br/></div><div class='line' id='LC279'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">antimuon_dres</span><span class="p">.</span><span class="n">lepton_photon</span><span class="p">.</span><span class="n">empty</span><span class="p">()</span> <span class="o">&amp;&amp;</span> <span class="n">muon_p</span><span class="p">.</span><span class="n">pT</span><span class="p">()</span> <span class="o">&gt;</span> <span class="mi">20</span><span class="p">)</span></div><div class='line' id='LC280'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC281'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">antimuon_dres</span><span class="p">.</span><span class="n">p_part</span> <span class="o">=</span> <span class="n">muon_p</span><span class="p">;</span></div><div class='line' id='LC282'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">antimuon_dres</span><span class="p">.</span><span class="n">lepton_photon</span> <span class="o">=</span> <span class="n">antimuon_photons</span><span class="p">;</span></div><div class='line' id='LC283'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC284'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">else</span></div><div class='line' id='LC285'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC286'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">muon_p</span><span class="p">.</span><span class="n">pT</span><span class="p">()</span> <span class="o">&gt;</span> <span class="n">antimuon_dres</span><span class="p">.</span><span class="n">p_part</span><span class="p">.</span><span class="n">pT</span><span class="p">()</span></div><div class='line' id='LC287'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="o">&amp;&amp;</span> <span class="n">muon_p</span><span class="p">.</span><span class="n">pT</span><span class="p">()</span> <span class="o">&gt;</span> <span class="mi">20</span><span class="p">)</span></div><div class='line' id='LC288'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC289'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">antimuon_dres</span><span class="p">.</span><span class="n">p_part</span> <span class="o">=</span> <span class="n">muon_p</span><span class="p">;</span></div><div class='line' id='LC290'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">antimuon_dres</span><span class="p">.</span><span class="n">lepton_photon</span> <span class="o">=</span> <span class="n">antimuon_photons</span><span class="p">;</span></div><div class='line' id='LC291'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC292'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC293'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC294'><br/></div><div class='line' id='LC295'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC296'><br/></div><div class='line' id='LC297'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// construct the Z candidates momentum</span></div><div class='line' id='LC298'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">FourMomentum</span> <span class="n">Z_momentum_ele</span><span class="p">(</span><span class="n">add</span><span class="p">(</span><span class="n">ele_dres</span><span class="p">.</span><span class="n">p_part</span><span class="p">,</span> <span class="n">pos_dres</span><span class="p">.</span><span class="n">p_part</span><span class="p">));</span></div><div class='line' id='LC299'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">FourMomentum</span> <span class="nf">Z_momentum_muon</span><span class="p">(</span><span class="n">add</span><span class="p">(</span><span class="n">muon_dres</span><span class="p">.</span><span class="n">p_part</span><span class="p">,</span> <span class="n">antimuon_dres</span><span class="p">.</span><span class="n">p_part</span><span class="p">));</span></div><div class='line' id='LC300'><br/></div><div class='line' id='LC301'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span><span class="p">(</span><span class="n">Z_momentum_ele</span><span class="p">.</span><span class="n">mass2</span><span class="p">()</span><span class="o">&lt;</span><span class="mi">0</span><span class="p">){</span></div><div class='line' id='LC302'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">MSG_WARNING</span><span class="p">(</span><span class="s">&quot;WARNING: negative e+e- invariant mass squared&quot;</span> <span class="o">&lt;&lt;</span> <span class="n">Z_momentum_ele</span><span class="p">.</span><span class="n">mass2</span><span class="p">()</span> <span class="p">);</span></div><div class='line' id='LC303'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">vetoEvent</span><span class="p">;</span></div><div class='line' id='LC304'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC305'><br/></div><div class='line' id='LC306'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span><span class="p">(</span><span class="n">Z_momentum_muon</span><span class="p">.</span><span class="n">mass2</span><span class="p">()</span><span class="o">&lt;</span><span class="mi">0</span><span class="p">){</span></div><div class='line' id='LC307'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">MSG_WARNING</span><span class="p">(</span><span class="s">&quot;WARNING: negative mu+mu- invariant mass squared&quot;</span> <span class="o">&lt;&lt;</span> <span class="n">Z_momentum_muon</span><span class="p">.</span><span class="n">mass2</span><span class="p">());</span></div><div class='line' id='LC308'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">vetoEvent</span><span class="p">;</span></div><div class='line' id='LC309'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC310'><br/></div><div class='line' id='LC311'><br/></div><div class='line' id='LC312'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">l1_eta</span><span class="p">(</span><span class="o">-</span><span class="mf">999.</span><span class="p">),</span><span class="n">l1_phi</span><span class="p">(</span><span class="o">-</span><span class="mf">999.</span><span class="p">),</span><span class="n">l2_eta</span><span class="p">(</span><span class="o">-</span><span class="mf">999.</span><span class="p">),</span><span class="n">l2_phi</span><span class="p">(</span><span class="o">-</span><span class="mf">999.</span><span class="p">);</span></div><div class='line' id='LC313'><br/></div><div class='line' id='LC314'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// Z mass window</span></div><div class='line' id='LC315'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">((</span><span class="n">Z_momentum_ele</span><span class="p">.</span><span class="n">mass</span><span class="p">()</span> <span class="o">&gt;</span> <span class="mf">71.</span> <span class="o">&amp;&amp;</span> <span class="n">Z_momentum_ele</span><span class="p">.</span><span class="n">mass</span><span class="p">()</span> <span class="o">&lt;</span> <span class="mf">111.</span><span class="p">)</span> <span class="o">&amp;&amp;</span> <span class="p">(</span><span class="o">!</span><span class="n">ele_dres</span><span class="p">.</span><span class="n">lepton_photon</span><span class="p">.</span><span class="n">empty</span><span class="p">()</span> <span class="o">&amp;&amp;</span> <span class="o">!</span><span class="n">pos_dres</span><span class="p">.</span><span class="n">lepton_photon</span><span class="p">.</span><span class="n">empty</span><span class="p">()))</span></div><div class='line' id='LC316'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC317'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">Z_ele</span> <span class="o">=</span> <span class="nb">true</span><span class="p">;</span></div><div class='line' id='LC318'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">l1_eta</span> <span class="o">=</span> <span class="n">ele_dres</span><span class="p">.</span><span class="n">p_part</span><span class="p">.</span><span class="n">eta</span><span class="p">();</span></div><div class='line' id='LC319'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">l1_phi</span> <span class="o">=</span> <span class="n">ele_dres</span><span class="p">.</span><span class="n">p_part</span><span class="p">.</span><span class="n">phi</span><span class="p">();</span></div><div class='line' id='LC320'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">l2_eta</span> <span class="o">=</span> <span class="n">pos_dres</span><span class="p">.</span><span class="n">p_part</span><span class="p">.</span><span class="n">eta</span><span class="p">();</span></div><div class='line' id='LC321'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">l2_phi</span> <span class="o">=</span> <span class="n">pos_dres</span><span class="p">.</span><span class="n">p_part</span><span class="p">.</span><span class="n">phi</span><span class="p">();</span></div><div class='line' id='LC322'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC323'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">((</span><span class="n">Z_momentum_muon</span><span class="p">.</span><span class="n">mass</span><span class="p">()</span> <span class="o">&gt;</span> <span class="mf">71.</span> <span class="o">&amp;&amp;</span> <span class="n">Z_momentum_muon</span><span class="p">.</span><span class="n">mass</span><span class="p">()</span> <span class="o">&lt;</span> <span class="mf">111.</span><span class="p">)</span> <span class="o">&amp;&amp;</span> <span class="p">(</span><span class="o">!</span><span class="n">muon_dres</span><span class="p">.</span><span class="n">lepton_photon</span><span class="p">.</span><span class="n">empty</span><span class="p">()</span> <span class="o">&amp;&amp;</span> <span class="o">!</span><span class="n">antimuon_dres</span><span class="p">.</span><span class="n">lepton_photon</span><span class="p">.</span><span class="n">empty</span><span class="p">()))</span></div><div class='line' id='LC324'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC325'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">Z_muon</span> <span class="o">=</span> <span class="nb">true</span><span class="p">;</span></div><div class='line' id='LC326'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">l1_eta</span> <span class="o">=</span> <span class="n">muon_dres</span><span class="p">.</span><span class="n">p_part</span><span class="p">.</span><span class="n">eta</span><span class="p">();</span></div><div class='line' id='LC327'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">l1_phi</span> <span class="o">=</span> <span class="n">muon_dres</span><span class="p">.</span><span class="n">p_part</span><span class="p">.</span><span class="n">phi</span><span class="p">();</span></div><div class='line' id='LC328'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">l2_eta</span> <span class="o">=</span> <span class="n">antimuon_dres</span><span class="p">.</span><span class="n">p_part</span><span class="p">.</span><span class="n">eta</span><span class="p">();</span></div><div class='line' id='LC329'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">l2_phi</span> <span class="o">=</span> <span class="n">antimuon_dres</span><span class="p">.</span><span class="n">p_part</span><span class="p">.</span><span class="n">phi</span><span class="p">();</span></div><div class='line' id='LC330'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC331'><br/></div><div class='line' id='LC332'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="o">!</span><span class="n">Z_ele</span> <span class="o">&amp;&amp;</span> <span class="o">!</span><span class="n">Z_muon</span><span class="p">)</span></div><div class='line' id='LC333'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">vetoEvent</span><span class="p">;</span></div><div class='line' id='LC334'><br/></div><div class='line' id='LC335'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">Z_ele</span> <span class="o">&amp;&amp;</span> <span class="n">Z_muon</span><span class="p">)</span></div><div class='line' id='LC336'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">vetoEvent</span><span class="p">;</span></div><div class='line' id='LC337'><br/></div><div class='line' id='LC338'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// particles to be removed from jets</span></div><div class='line' id='LC339'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">vector</span><span class="o">&lt;</span><span class="kt">unsigned</span> <span class="kt">int</span><span class="o">&gt;</span> <span class="n">canc_part</span><span class="p">;</span></div><div class='line' id='LC340'><br/></div><div class='line' id='LC341'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">Z_muon</span><span class="p">)</span></div><div class='line' id='LC342'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC343'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">canc_part</span> <span class="o">=</span> <span class="n">muon_dres</span><span class="p">.</span><span class="n">lepton_photon</span><span class="p">;</span></div><div class='line' id='LC344'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">for</span> <span class="p">(</span><span class="kt">unsigned</span> <span class="kt">int</span> <span class="n">l</span> <span class="o">=</span> <span class="mi">0</span><span class="p">;</span> <span class="n">l</span> <span class="o">&lt;</span> <span class="n">antimuon_dres</span><span class="p">.</span><span class="n">lepton_photon</span><span class="p">.</span><span class="n">size</span><span class="p">();</span> <span class="n">l</span><span class="o">++</span><span class="p">)</span></div><div class='line' id='LC345'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">canc_part</span><span class="p">.</span><span class="n">push_back</span><span class="p">(</span><span class="n">antimuon_dres</span><span class="p">.</span><span class="n">lepton_photon</span><span class="p">.</span><span class="n">at</span><span class="p">(</span><span class="n">l</span><span class="p">));</span></div><div class='line' id='LC346'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC347'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">else</span></div><div class='line' id='LC348'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC349'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">canc_part</span> <span class="o">=</span> <span class="n">ele_dres</span><span class="p">.</span><span class="n">lepton_photon</span><span class="p">;</span></div><div class='line' id='LC350'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">for</span> <span class="p">(</span><span class="kt">unsigned</span> <span class="kt">int</span> <span class="n">y</span> <span class="o">=</span> <span class="mi">0</span><span class="p">;</span> <span class="n">y</span> <span class="o">&lt;</span> <span class="n">pos_dres</span><span class="p">.</span><span class="n">lepton_photon</span><span class="p">.</span><span class="n">size</span><span class="p">();</span> <span class="n">y</span><span class="o">++</span><span class="p">)</span></div><div class='line' id='LC351'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">canc_part</span><span class="p">.</span><span class="n">push_back</span><span class="p">(</span><span class="n">pos_dres</span><span class="p">.</span><span class="n">lepton_photon</span><span class="p">.</span><span class="n">at</span><span class="p">(</span><span class="n">y</span><span class="p">));</span></div><div class='line' id='LC352'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC353'><br/></div><div class='line' id='LC354'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="o">!</span><span class="n">neutrinos</span><span class="p">.</span><span class="n">empty</span><span class="p">())</span></div><div class='line' id='LC355'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC356'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">for</span> <span class="p">(</span><span class="kt">unsigned</span> <span class="kt">int</span> <span class="n">k</span> <span class="o">=</span> <span class="mi">0</span><span class="p">;</span> <span class="n">k</span> <span class="o">&lt;</span> <span class="n">neutrinos</span><span class="p">.</span><span class="n">size</span><span class="p">();</span> <span class="n">k</span><span class="o">++</span><span class="p">)</span></div><div class='line' id='LC357'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">canc_part</span><span class="p">.</span><span class="n">push_back</span><span class="p">(</span><span class="n">neutrinos</span><span class="p">.</span><span class="n">at</span><span class="p">(</span><span class="n">k</span><span class="p">));</span></div><div class='line' id='LC358'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC359'><br/></div><div class='line' id='LC360'><br/></div><div class='line' id='LC361'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">stable_sort</span><span class="p">(</span><span class="n">canc_part</span><span class="p">.</span><span class="n">begin</span><span class="p">(),</span> <span class="n">canc_part</span><span class="p">.</span><span class="n">end</span><span class="p">());</span></div><div class='line' id='LC362'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">for</span> <span class="p">(</span><span class="kt">int</span> <span class="n">i</span> <span class="o">=</span> <span class="p">(</span><span class="n">canc_part</span><span class="p">.</span><span class="n">size</span><span class="p">()</span> <span class="o">-</span> <span class="mi">1</span><span class="p">);</span> <span class="n">i</span> <span class="o">&gt;=</span> <span class="mi">0</span><span class="p">;</span> <span class="n">i</span><span class="o">--</span><span class="p">)</span></div><div class='line' id='LC363'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">part_jets</span><span class="p">.</span><span class="n">erase</span><span class="p">(</span><span class="n">part_jets</span><span class="p">.</span><span class="n">begin</span><span class="p">()</span> <span class="o">+</span> <span class="n">canc_part</span><span class="p">.</span><span class="n">at</span><span class="p">(</span><span class="n">i</span><span class="p">));</span></div><div class='line' id='LC364'><br/></div><div class='line' id='LC365'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// jet clustering</span></div><div class='line' id='LC366'><br/></div><div class='line' id='LC367'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">for</span> <span class="p">(</span><span class="kt">unsigned</span> <span class="kt">int</span> <span class="n">l</span> <span class="o">=</span> <span class="mi">0</span><span class="p">;</span> <span class="n">l</span> <span class="o">&lt;</span> <span class="n">part_jets</span><span class="p">.</span><span class="n">size</span><span class="p">();</span> <span class="n">l</span><span class="o">++</span><span class="p">)</span></div><div class='line' id='LC368'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC369'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">FourMomentum</span> <span class="n">p</span><span class="p">(</span><span class="n">part_jets</span><span class="p">.</span><span class="n">at</span><span class="p">(</span><span class="n">l</span><span class="p">).</span><span class="n">momentum</span><span class="p">());</span></div><div class='line' id='LC370'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">fastjet</span><span class="o">::</span><span class="n">PseudoJet</span> <span class="n">pseudoJet</span><span class="p">(</span><span class="n">p</span><span class="p">.</span><span class="n">px</span><span class="p">(),</span> <span class="n">p</span><span class="p">.</span><span class="n">py</span><span class="p">(),</span> <span class="n">p</span><span class="p">.</span><span class="n">pz</span><span class="p">(),</span> <span class="n">p</span><span class="p">.</span><span class="n">E</span><span class="p">());</span></div><div class='line' id='LC371'><br/></div><div class='line' id='LC372'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">pseudoJet</span><span class="p">.</span><span class="n">set_user_index</span><span class="p">(</span><span class="n">l</span><span class="p">);</span></div><div class='line' id='LC373'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">vecs</span><span class="p">.</span><span class="n">push_back</span><span class="p">(</span><span class="n">pseudoJet</span><span class="p">);</span></div><div class='line' id='LC374'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC375'><br/></div><div class='line' id='LC376'><br/></div><div class='line' id='LC377'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// jet selection</span></div><div class='line' id='LC378'><br/></div><div class='line' id='LC379'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">ht</span> <span class="o">=</span> <span class="mi">0</span><span class="p">;</span></div><div class='line' id='LC380'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">vector</span><span class="o">&lt;</span><span class="n">fastjet</span><span class="o">::</span><span class="n">PseudoJet</span><span class="o">&gt;</span> <span class="n">jet_list</span><span class="p">;</span></div><div class='line' id='LC381'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">fastjet</span><span class="o">::</span><span class="n">ClusterSequence</span> <span class="n">cseq</span><span class="p">(</span><span class="n">vecs</span><span class="p">,</span> <span class="n">fastjet</span><span class="o">::</span><span class="n">JetDefinition</span><span class="p">(</span><span class="n">fastjet</span><span class="o">::</span> <span class="n">antikt_algorithm</span><span class="p">,</span> <span class="mf">0.5</span><span class="p">));</span></div><div class='line' id='LC382'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">vector</span><span class="o">&lt;</span><span class="n">fastjet</span><span class="o">::</span><span class="n">PseudoJet</span><span class="o">&gt;</span> <span class="n">jets</span> <span class="o">=</span> <span class="n">sorted_by_pt</span><span class="p">(</span><span class="n">cseq</span><span class="p">.</span><span class="n">inclusive_jets</span><span class="p">(</span><span class="mf">30.0</span><span class="p">));</span></div><div class='line' id='LC383'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">for</span> <span class="p">(</span><span class="kt">unsigned</span> <span class="kt">int</span> <span class="n">i</span> <span class="o">=</span> <span class="mi">0</span><span class="p">;</span> <span class="n">i</span> <span class="o">&lt;</span> <span class="n">jets</span><span class="p">.</span><span class="n">size</span><span class="p">();</span> <span class="n">i</span><span class="o">++</span><span class="p">)</span></div><div class='line' id='LC384'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC385'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">etaj</span> <span class="o">=</span> <span class="n">jets</span><span class="p">[</span><span class="n">i</span><span class="p">].</span><span class="n">eta</span><span class="p">();</span></div><div class='line' id='LC386'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">phij</span> <span class="o">=</span> <span class="n">jets</span><span class="p">[</span><span class="n">i</span><span class="p">].</span><span class="n">phi</span><span class="p">();</span></div><div class='line' id='LC387'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">bool</span> <span class="n">found</span> <span class="o">=</span> <span class="nb">false</span><span class="p">;</span></div><div class='line' id='LC388'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">vector</span><span class="o">&lt;</span><span class="n">fastjet</span><span class="o">::</span><span class="n">PseudoJet</span><span class="o">&gt;</span> <span class="n">constituents</span> <span class="o">=</span> <span class="n">cseq</span><span class="p">.</span><span class="n">constituents</span><span class="p">(</span><span class="n">jets</span><span class="p">[</span><span class="n">i</span><span class="p">]);</span></div><div class='line' id='LC389'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">for</span> <span class="p">(</span><span class="kt">unsigned</span> <span class="kt">int</span> <span class="n">s</span> <span class="o">=</span> <span class="mi">0</span><span class="p">;</span> <span class="n">s</span> <span class="o">&lt;</span> <span class="n">constituents</span><span class="p">.</span><span class="n">size</span><span class="p">()</span> <span class="o">&amp;&amp;</span> <span class="o">!</span><span class="n">found</span><span class="p">;</span> <span class="n">s</span><span class="o">++</span><span class="p">)</span></div><div class='line' id='LC390'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC391'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">int</span> <span class="n">index</span> <span class="o">=</span> <span class="n">constituents</span><span class="p">.</span><span class="n">at</span><span class="p">(</span><span class="n">s</span><span class="p">).</span><span class="n">user_index</span><span class="p">();</span></div><div class='line' id='LC392'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">PID</span><span class="o">::</span><span class="n">threeCharge</span><span class="p">(</span><span class="n">part_jets</span><span class="p">.</span><span class="n">at</span><span class="p">(</span><span class="n">index</span><span class="p">).</span><span class="n">pdgId</span><span class="p">())</span> <span class="o">!=</span> <span class="mi">0</span><span class="p">);</span></div><div class='line' id='LC393'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">found</span> <span class="o">=</span> <span class="nb">true</span><span class="p">;</span></div><div class='line' id='LC394'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC395'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">fabs</span><span class="p">(</span><span class="n">etaj</span><span class="p">)</span> <span class="o">&lt;</span> <span class="mf">2.4</span> <span class="o">&amp;&amp;</span> <span class="n">found</span> <span class="o">&amp;&amp;</span> <span class="n">jets</span><span class="p">[</span><span class="n">i</span><span class="p">].</span><span class="n">perp</span><span class="p">()</span> <span class="o">&gt;</span> <span class="mi">30</span><span class="p">)</span></div><div class='line' id='LC396'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC397'><br/></div><div class='line' id='LC398'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">DR1</span> <span class="o">=</span> <span class="n">Rivet</span><span class="o">::</span><span class="n">deltaR</span><span class="p">(</span><span class="n">etaj</span><span class="p">,</span><span class="n">phij</span><span class="p">,</span><span class="n">l1_eta</span><span class="p">,</span><span class="n">l1_phi</span><span class="p">);</span></div><div class='line' id='LC399'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">DR2</span> <span class="o">=</span> <span class="n">Rivet</span><span class="o">::</span><span class="n">deltaR</span><span class="p">(</span><span class="n">etaj</span><span class="p">,</span><span class="n">phij</span><span class="p">,</span><span class="n">l2_eta</span><span class="p">,</span><span class="n">l2_phi</span><span class="p">);</span></div><div class='line' id='LC400'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span> <span class="n">std</span><span class="o">::</span><span class="n">min</span><span class="p">(</span><span class="n">DR1</span><span class="p">,</span><span class="n">DR2</span><span class="p">)</span> <span class="o">&lt;</span> <span class="mf">0.5</span> <span class="p">)</span> <span class="p">{</span> </div><div class='line' id='LC401'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// drop jet closer than 0.5 from the selected collection</span></div><div class='line' id='LC402'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">continue</span><span class="p">;</span></div><div class='line' id='LC403'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// drop event if at least one jet close than 0.5 in the selected collection</span></div><div class='line' id='LC404'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">//                vetoEvent;</span></div><div class='line' id='LC405'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC406'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC407'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">jet_list</span><span class="p">.</span><span class="n">push_back</span><span class="p">(</span><span class="n">jets</span><span class="p">[</span><span class="n">i</span><span class="p">]);</span></div><div class='line' id='LC408'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">ht</span> <span class="o">=</span> <span class="n">ht</span> <span class="o">+</span> <span class="n">jets</span><span class="p">[</span><span class="n">i</span><span class="p">].</span><span class="n">perp</span><span class="p">();</span></div><div class='line' id='LC409'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC410'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC411'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC412'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">jet_list</span><span class="p">.</span><span class="n">size</span><span class="p">()</span> <span class="o">&gt;</span> <span class="mi">0</span><span class="p">)</span></div><div class='line' id='LC413'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC414'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">leadingjet_pt</span> <span class="o">=</span> <span class="n">jets</span><span class="p">[</span><span class="mi">0</span><span class="p">].</span><span class="n">perp</span><span class="p">();</span></div><div class='line' id='LC415'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">leadingEta</span> <span class="o">=</span> <span class="n">jets</span><span class="p">[</span><span class="mi">0</span><span class="p">].</span><span class="n">eta</span><span class="p">();</span></div><div class='line' id='LC416'><br/></div><div class='line' id='LC417'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_leading_jet_pt_tot</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">leadingjet_pt</span><span class="p">,</span> <span class="n">weight</span><span class="o">*</span><span class="mf">0.5</span><span class="p">);</span></div><div class='line' id='LC418'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_leading_jet_eta_tot</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">std</span><span class="o">::</span><span class="n">fabs</span><span class="p">(</span><span class="n">leadingEta</span><span class="p">),</span> <span class="n">weight</span><span class="o">*</span><span class="mf">0.5</span><span class="p">);</span></div><div class='line' id='LC419'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_ht1_tot</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">ht</span><span class="p">,</span> <span class="n">weight</span><span class="o">*</span><span class="mf">0.5</span><span class="p">);</span></div><div class='line' id='LC420'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_excmult_jets_tot</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">jet_list</span><span class="p">.</span><span class="n">size</span><span class="p">(),</span> <span class="n">weight</span><span class="o">*</span><span class="mf">0.5</span><span class="p">);</span></div><div class='line' id='LC421'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">for</span> <span class="p">(</span> <span class="kt">unsigned</span> <span class="kt">int</span> <span class="n">iter</span> <span class="o">=</span> <span class="mi">1</span><span class="p">;</span> <span class="n">iter</span> <span class="o">&lt;=</span> <span class="n">jet_list</span><span class="p">.</span><span class="n">size</span><span class="p">();</span> <span class="n">iter</span><span class="o">++</span> <span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC422'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_incmult_jets_tot</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">iter</span><span class="p">,</span> <span class="n">weight</span><span class="o">*</span><span class="mf">0.5</span><span class="p">);</span></div><div class='line' id='LC423'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC424'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">Z_ele</span><span class="p">)</span></div><div class='line' id='LC425'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC426'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_leading_jet_pt_ele</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">leadingjet_pt</span><span class="p">,</span> <span class="n">weight</span><span class="p">);</span></div><div class='line' id='LC427'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_leading_jet_eta_ele</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">std</span><span class="o">::</span><span class="n">fabs</span><span class="p">(</span><span class="n">leadingEta</span><span class="p">),</span> <span class="n">weight</span><span class="p">);</span></div><div class='line' id='LC428'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_excmult_jets_ele</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">jet_list</span><span class="p">.</span><span class="n">size</span><span class="p">(),</span> <span class="n">weight</span><span class="p">);</span></div><div class='line' id='LC429'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_ht1_ele</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">ht</span><span class="p">,</span> <span class="n">weight</span><span class="p">);</span></div><div class='line' id='LC430'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">for</span> <span class="p">(</span> <span class="kt">unsigned</span> <span class="kt">int</span> <span class="n">iter</span> <span class="o">=</span> <span class="mi">1</span><span class="p">;</span> <span class="n">iter</span> <span class="o">&lt;=</span> <span class="n">jet_list</span><span class="p">.</span><span class="n">size</span><span class="p">();</span> <span class="n">iter</span><span class="o">++</span> <span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC431'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_incmult_jets_ele</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">iter</span><span class="p">,</span> <span class="n">weight</span><span class="p">);</span></div><div class='line' id='LC432'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC433'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC434'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">Z_muon</span><span class="p">)</span></div><div class='line' id='LC435'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC436'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_leading_jet_pt_muo</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">leadingjet_pt</span><span class="p">,</span> <span class="n">weight</span><span class="p">);</span></div><div class='line' id='LC437'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_leading_jet_eta_muo</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">std</span><span class="o">::</span><span class="n">fabs</span><span class="p">(</span><span class="n">leadingEta</span><span class="p">),</span> <span class="n">weight</span><span class="p">);</span></div><div class='line' id='LC438'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_excmult_jets_muo</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">jet_list</span><span class="p">.</span><span class="n">size</span><span class="p">(),</span> <span class="n">weight</span><span class="p">);</span></div><div class='line' id='LC439'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_ht1_muo</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">ht</span><span class="p">,</span> <span class="n">weight</span><span class="p">);</span></div><div class='line' id='LC440'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">for</span> <span class="p">(</span> <span class="kt">unsigned</span> <span class="kt">int</span> <span class="n">iter</span> <span class="o">=</span> <span class="mi">1</span><span class="p">;</span> <span class="n">iter</span> <span class="o">&lt;=</span> <span class="n">jet_list</span><span class="p">.</span><span class="n">size</span><span class="p">();</span> <span class="n">iter</span><span class="o">++</span> <span class="p">)</span> <span class="p">{</span></div><div class='line' id='LC441'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_incmult_jets_muo</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">iter</span><span class="p">,</span> <span class="n">weight</span><span class="p">);</span></div><div class='line' id='LC442'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC443'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC444'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC445'><br/></div><div class='line' id='LC446'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">jet_list</span><span class="p">.</span><span class="n">size</span><span class="p">()</span> <span class="o">&gt;</span> <span class="mi">1</span><span class="p">)</span></div><div class='line' id='LC447'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC448'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">second_jet_pt</span> <span class="o">=</span> <span class="n">jets</span><span class="p">[</span><span class="mi">1</span><span class="p">].</span><span class="n">perp</span><span class="p">();</span></div><div class='line' id='LC449'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">second_jet_eta</span> <span class="o">=</span> <span class="n">jets</span><span class="p">[</span><span class="mi">1</span><span class="p">].</span><span class="n">eta</span><span class="p">();</span></div><div class='line' id='LC450'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_second_jet_pt_tot</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">second_jet_pt</span><span class="p">,</span> <span class="n">weight</span><span class="o">*</span><span class="mf">0.5</span><span class="p">);</span></div><div class='line' id='LC451'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_second_jet_eta_tot</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">std</span><span class="o">::</span><span class="n">fabs</span><span class="p">(</span><span class="n">second_jet_eta</span><span class="p">),</span> <span class="n">weight</span><span class="o">*</span><span class="mf">0.5</span><span class="p">);</span></div><div class='line' id='LC452'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_ht2_tot</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">ht</span><span class="p">,</span> <span class="n">weight</span><span class="o">*</span><span class="mf">0.5</span><span class="p">);</span></div><div class='line' id='LC453'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">Z_ele</span><span class="p">)</span></div><div class='line' id='LC454'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC455'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_second_jet_pt_ele</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">second_jet_pt</span><span class="p">,</span> <span class="n">weight</span><span class="p">);</span></div><div class='line' id='LC456'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_second_jet_eta_ele</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">std</span><span class="o">::</span><span class="n">fabs</span><span class="p">(</span><span class="n">second_jet_eta</span><span class="p">),</span> <span class="n">weight</span><span class="p">);</span></div><div class='line' id='LC457'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_ht2_ele</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">ht</span><span class="p">,</span> <span class="n">weight</span><span class="p">);</span></div><div class='line' id='LC458'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC459'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">Z_muon</span><span class="p">)</span></div><div class='line' id='LC460'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC461'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_second_jet_pt_muo</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">second_jet_pt</span><span class="p">,</span> <span class="n">weight</span><span class="p">);</span></div><div class='line' id='LC462'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_second_jet_eta_muo</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">std</span><span class="o">::</span><span class="n">fabs</span><span class="p">(</span><span class="n">second_jet_eta</span><span class="p">),</span> <span class="n">weight</span><span class="p">);</span></div><div class='line' id='LC463'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_ht2_muo</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">ht</span><span class="p">,</span> <span class="n">weight</span><span class="p">);</span></div><div class='line' id='LC464'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC465'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC466'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">jet_list</span><span class="p">.</span><span class="n">size</span><span class="p">()</span> <span class="o">&gt;</span> <span class="mi">2</span><span class="p">)</span></div><div class='line' id='LC467'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC468'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">third_jet_pt</span> <span class="o">=</span> <span class="n">jets</span><span class="p">[</span><span class="mi">2</span><span class="p">].</span><span class="n">perp</span><span class="p">();</span></div><div class='line' id='LC469'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">third_jet_eta</span> <span class="o">=</span> <span class="n">jets</span><span class="p">[</span><span class="mi">2</span><span class="p">].</span><span class="n">eta</span><span class="p">();</span></div><div class='line' id='LC470'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_third_jet_pt_tot</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">third_jet_pt</span><span class="p">,</span> <span class="n">weight</span><span class="o">*</span><span class="mf">0.5</span><span class="p">);</span></div><div class='line' id='LC471'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_third_jet_eta_tot</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">std</span><span class="o">::</span><span class="n">fabs</span><span class="p">(</span><span class="n">third_jet_eta</span><span class="p">),</span> <span class="n">weight</span><span class="o">*</span><span class="mf">0.5</span><span class="p">);</span></div><div class='line' id='LC472'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_ht3_tot</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">ht</span><span class="p">,</span> <span class="n">weight</span><span class="o">*</span><span class="mf">0.5</span><span class="p">);</span></div><div class='line' id='LC473'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">Z_ele</span><span class="p">)</span></div><div class='line' id='LC474'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC475'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_third_jet_pt_ele</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">third_jet_pt</span><span class="p">,</span> <span class="n">weight</span><span class="p">);</span></div><div class='line' id='LC476'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_third_jet_eta_ele</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">std</span><span class="o">::</span><span class="n">fabs</span><span class="p">(</span><span class="n">third_jet_eta</span><span class="p">),</span> <span class="n">weight</span><span class="p">);</span></div><div class='line' id='LC477'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_ht3_ele</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">ht</span><span class="p">,</span> <span class="n">weight</span><span class="p">);</span></div><div class='line' id='LC478'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC479'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">Z_muon</span><span class="p">)</span></div><div class='line' id='LC480'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC481'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_third_jet_pt_muo</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">third_jet_pt</span><span class="p">,</span> <span class="n">weight</span><span class="p">);</span></div><div class='line' id='LC482'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_third_jet_eta_muo</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">std</span><span class="o">::</span><span class="n">fabs</span><span class="p">(</span><span class="n">third_jet_eta</span><span class="p">),</span> <span class="n">weight</span><span class="p">);</span></div><div class='line' id='LC483'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_ht3_muo</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">ht</span><span class="p">,</span> <span class="n">weight</span><span class="p">);</span></div><div class='line' id='LC484'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC485'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC486'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">jet_list</span><span class="p">.</span><span class="n">size</span><span class="p">()</span> <span class="o">&gt;</span> <span class="mi">3</span><span class="p">)</span></div><div class='line' id='LC487'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC488'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">fourth_jet_pt</span> <span class="o">=</span> <span class="n">jets</span><span class="p">[</span><span class="mi">3</span><span class="p">].</span><span class="n">perp</span><span class="p">();</span></div><div class='line' id='LC489'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">fourth_jet_eta</span> <span class="o">=</span> <span class="n">jets</span><span class="p">[</span><span class="mi">3</span><span class="p">].</span><span class="n">eta</span><span class="p">();</span></div><div class='line' id='LC490'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_fourth_jet_pt_tot</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">fourth_jet_pt</span><span class="p">,</span> <span class="n">weight</span><span class="o">*</span><span class="mf">0.5</span><span class="p">);</span></div><div class='line' id='LC491'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_fourth_jet_eta_tot</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">std</span><span class="o">::</span><span class="n">fabs</span><span class="p">(</span><span class="n">fourth_jet_eta</span><span class="p">),</span> <span class="n">weight</span><span class="o">*</span><span class="mf">0.5</span><span class="p">);</span></div><div class='line' id='LC492'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_ht4_tot</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">ht</span><span class="p">,</span> <span class="n">weight</span><span class="o">*</span><span class="mf">0.5</span><span class="p">);</span></div><div class='line' id='LC493'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">Z_ele</span><span class="p">)</span></div><div class='line' id='LC494'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC495'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_fourth_jet_pt_ele</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">fourth_jet_pt</span><span class="p">,</span> <span class="n">weight</span><span class="p">);</span></div><div class='line' id='LC496'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_fourth_jet_eta_ele</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">std</span><span class="o">::</span><span class="n">fabs</span><span class="p">(</span><span class="n">fourth_jet_eta</span><span class="p">),</span> <span class="n">weight</span><span class="p">);</span></div><div class='line' id='LC497'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_ht4_ele</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">ht</span><span class="p">,</span> <span class="n">weight</span><span class="p">);</span></div><div class='line' id='LC498'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC499'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">Z_muon</span><span class="p">)</span></div><div class='line' id='LC500'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC501'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_fourth_jet_pt_muo</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">fourth_jet_pt</span><span class="p">,</span> <span class="n">weight</span><span class="p">);</span></div><div class='line' id='LC502'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_fourth_jet_eta_muo</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">std</span><span class="o">::</span><span class="n">fabs</span><span class="p">(</span><span class="n">fourth_jet_eta</span><span class="p">),</span> <span class="n">weight</span><span class="p">);</span></div><div class='line' id='LC503'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">_h_ht4_muo</span><span class="o">-&gt;</span><span class="n">fill</span><span class="p">(</span><span class="n">ht</span><span class="p">,</span> <span class="n">weight</span><span class="p">);</span></div><div class='line' id='LC504'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC505'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC506'><br/></div><div class='line' id='LC507'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC508'><br/></div><div class='line' id='LC509'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">/// Normalise histograms etc., after the run</span></div><div class='line' id='LC510'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">void</span> <span class="n">finalize</span><span class="p">()</span></div><div class='line' id='LC511'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">{</span></div><div class='line' id='LC512'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC513'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">MSG_INFO</span><span class="p">(</span><span class="s">&quot;Cross section = &quot;</span> <span class="o">&lt;&lt;</span> <span class="n">std</span><span class="o">::</span><span class="n">setfill</span><span class="p">(</span><span class="sc">&#39; &#39;</span><span class="p">)</span> <span class="o">&lt;&lt;</span> <span class="n">std</span><span class="o">::</span><span class="n">setw</span><span class="p">(</span><span class="mi">14</span><span class="p">)</span> <span class="o">&lt;&lt;</span> <span class="n">std</span><span class="o">::</span><span class="n">fixed</span> <span class="o">&lt;&lt;</span> <span class="n">std</span><span class="o">::</span><span class="n">setprecision</span><span class="p">(</span><span class="mi">3</span><span class="p">)</span> <span class="o">&lt;&lt;</span> <span class="n">crossSection</span><span class="p">()</span> <span class="o">&lt;&lt;</span> <span class="s">&quot; pb&quot;</span><span class="p">);</span></div><div class='line' id='LC514'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">MSG_INFO</span><span class="p">(</span><span class="s">&quot;# Events      = &quot;</span> <span class="o">&lt;&lt;</span> <span class="n">std</span><span class="o">::</span><span class="n">setfill</span><span class="p">(</span><span class="sc">&#39; &#39;</span><span class="p">)</span> <span class="o">&lt;&lt;</span> <span class="n">std</span><span class="o">::</span><span class="n">setw</span><span class="p">(</span><span class="mi">14</span><span class="p">)</span> <span class="o">&lt;&lt;</span> <span class="n">std</span><span class="o">::</span><span class="n">fixed</span> <span class="o">&lt;&lt;</span> <span class="n">std</span><span class="o">::</span><span class="n">setprecision</span><span class="p">(</span><span class="mi">3</span><span class="p">)</span> <span class="o">&lt;&lt;</span> <span class="n">numEvents</span><span class="p">()</span> <span class="p">);</span></div><div class='line' id='LC515'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">MSG_INFO</span><span class="p">(</span><span class="s">&quot;SumW          = &quot;</span> <span class="o">&lt;&lt;</span> <span class="n">std</span><span class="o">::</span><span class="n">setfill</span><span class="p">(</span><span class="sc">&#39; &#39;</span><span class="p">)</span> <span class="o">&lt;&lt;</span> <span class="n">std</span><span class="o">::</span><span class="n">setw</span><span class="p">(</span><span class="mi">14</span><span class="p">)</span> <span class="o">&lt;&lt;</span> <span class="n">std</span><span class="o">::</span><span class="n">fixed</span> <span class="o">&lt;&lt;</span> <span class="n">std</span><span class="o">::</span><span class="n">setprecision</span><span class="p">(</span><span class="mi">3</span><span class="p">)</span> <span class="o">&lt;&lt;</span> <span class="n">sumOfWeights</span><span class="p">());</span></div><div class='line' id='LC516'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">MSG_INFO</span><span class="p">(</span><span class="s">&quot;SumW2         = &quot;</span> <span class="o">&lt;&lt;</span> <span class="n">std</span><span class="o">::</span><span class="n">setfill</span><span class="p">(</span><span class="sc">&#39; &#39;</span><span class="p">)</span> <span class="o">&lt;&lt;</span> <span class="n">std</span><span class="o">::</span><span class="n">setw</span><span class="p">(</span><span class="mi">14</span><span class="p">)</span> <span class="o">&lt;&lt;</span> <span class="n">std</span><span class="o">::</span><span class="n">fixed</span> <span class="o">&lt;&lt;</span> <span class="n">std</span><span class="o">::</span><span class="n">setprecision</span><span class="p">(</span><span class="mi">3</span><span class="p">)</span> <span class="o">&lt;&lt;</span> <span class="n">SumW2</span><span class="p">);</span></div><div class='line' id='LC517'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="nf">mean_weight</span><span class="p">(</span><span class="mf">0.</span><span class="p">);</span></div><div class='line' id='LC518'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span><span class="n">numEvents</span><span class="p">()</span> <span class="o">!=</span> <span class="mi">0</span><span class="p">)</span> <span class="n">mean_weight</span> <span class="o">=</span> <span class="n">sumOfWeights</span><span class="p">()</span><span class="o">/</span><span class="p">(</span><span class="kt">double</span><span class="p">)</span> <span class="n">numEvents</span><span class="p">();</span></div><div class='line' id='LC519'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">MSG_INFO</span><span class="p">(</span><span class="s">&quot;Mean weight   = &quot;</span> <span class="o">&lt;&lt;</span> <span class="n">std</span><span class="o">::</span><span class="n">setfill</span><span class="p">(</span><span class="sc">&#39; &#39;</span><span class="p">)</span> <span class="o">&lt;&lt;</span> <span class="n">std</span><span class="o">::</span><span class="n">setw</span><span class="p">(</span><span class="mi">14</span><span class="p">)</span> <span class="o">&lt;&lt;</span> <span class="n">std</span><span class="o">::</span><span class="n">fixed</span> <span class="o">&lt;&lt;</span> <span class="n">std</span><span class="o">::</span><span class="n">setprecision</span><span class="p">(</span><span class="mi">6</span><span class="p">)</span> <span class="o">&lt;&lt;</span> <span class="n">mean_weight</span><span class="p">);</span></div><div class='line' id='LC520'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="nf">norm</span><span class="p">(</span><span class="mf">1.</span><span class="p">);</span></div><div class='line' id='LC521'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="k">if</span> <span class="p">(</span> <span class="n">std</span><span class="o">::</span><span class="n">fabs</span><span class="p">(</span><span class="n">sumOfWeights</span><span class="p">())</span> <span class="o">&gt;</span> <span class="mf">0.</span> <span class="p">)</span> <span class="n">norm</span> <span class="o">=</span> <span class="n">crossSection</span><span class="p">()</span><span class="o">/</span><span class="n">sumOfWeights</span><span class="p">();</span></div><div class='line' id='LC522'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">MSG_INFO</span><span class="p">(</span><span class="s">&quot;Norm factor   = &quot;</span> <span class="o">&lt;&lt;</span> <span class="n">std</span><span class="o">::</span><span class="n">setfill</span><span class="p">(</span><span class="sc">&#39; &#39;</span><span class="p">)</span> <span class="o">&lt;&lt;</span> <span class="n">std</span><span class="o">::</span><span class="n">setw</span><span class="p">(</span><span class="mi">14</span><span class="p">)</span> <span class="o">&lt;&lt;</span> <span class="n">std</span><span class="o">::</span><span class="n">fixed</span> <span class="o">&lt;&lt;</span> <span class="n">std</span><span class="o">::</span><span class="n">setprecision</span><span class="p">(</span><span class="mi">6</span><span class="p">)</span> <span class="o">&lt;&lt;</span> <span class="n">norm</span><span class="p">);</span></div><div class='line' id='LC523'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</div><div class='line' id='LC524'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_excmult_jets_ele</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC525'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_incmult_jets_ele</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC526'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_leading_jet_pt_ele</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC527'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_second_jet_pt_ele</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC528'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_third_jet_pt_ele</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC529'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_fourth_jet_pt_ele</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC530'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_leading_jet_eta_ele</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC531'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_second_jet_eta_ele</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC532'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_third_jet_eta_ele</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC533'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_fourth_jet_eta_ele</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC534'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_ht1_ele</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC535'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_ht2_ele</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC536'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_ht3_ele</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC537'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_ht4_ele</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC538'><br/></div><div class='line' id='LC539'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_excmult_jets_muo</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC540'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_incmult_jets_muo</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC541'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_leading_jet_pt_muo</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC542'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_second_jet_pt_muo</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC543'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_third_jet_pt_muo</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC544'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_fourth_jet_pt_muo</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC545'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_leading_jet_eta_muo</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC546'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_second_jet_eta_muo</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC547'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_third_jet_eta_muo</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC548'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_fourth_jet_eta_muo</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC549'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_ht1_muo</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC550'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_ht2_muo</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC551'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_ht3_muo</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC552'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_ht4_muo</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC553'><br/></div><div class='line' id='LC554'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_excmult_jets_tot</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC555'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_incmult_jets_tot</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC556'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_leading_jet_pt_tot</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC557'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_second_jet_pt_tot</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC558'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_third_jet_pt_tot</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC559'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_fourth_jet_pt_tot</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC560'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_leading_jet_eta_tot</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC561'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_second_jet_eta_tot</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC562'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_third_jet_eta_tot</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC563'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_fourth_jet_eta_tot</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC564'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_ht1_tot</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC565'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_ht2_tot</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC566'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_ht3_tot</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC567'>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">scale</span><span class="p">(</span><span class="n">_h_ht4_tot</span><span class="p">,</span> <span class="n">norm</span> <span class="p">);</span></div><div class='line' id='LC568'><br/></div><div class='line' id='LC569'><br/></div><div class='line' id='LC570'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="p">}</span></div><div class='line' id='LC571'><br/></div><div class='line' id='LC572'>&nbsp;&nbsp;<span class="k">private</span><span class="o">:</span></div><div class='line' id='LC573'><br/></div><div class='line' id='LC574'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">// Data members like post-cuts event weight counters go here</span></div><div class='line' id='LC575'><br/></div><div class='line' id='LC576'>&nbsp;&nbsp;<span class="k">private</span><span class="o">:</span></div><div class='line' id='LC577'><br/></div><div class='line' id='LC578'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="kt">double</span> <span class="n">SumW2</span><span class="p">;</span></div><div class='line' id='LC579'><br/></div><div class='line' id='LC580'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c1">/// @name Histograms</span></div><div class='line' id='LC581'><br/></div><div class='line' id='LC582'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_excmult_jets_ele</span><span class="p">;</span></div><div class='line' id='LC583'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_incmult_jets_ele</span><span class="p">;</span></div><div class='line' id='LC584'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_leading_jet_pt_ele</span><span class="p">;</span></div><div class='line' id='LC585'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_second_jet_pt_ele</span><span class="p">;</span></div><div class='line' id='LC586'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_third_jet_pt_ele</span><span class="p">;</span></div><div class='line' id='LC587'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_fourth_jet_pt_ele</span><span class="p">;</span></div><div class='line' id='LC588'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_leading_jet_eta_ele</span><span class="p">;</span></div><div class='line' id='LC589'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_second_jet_eta_ele</span><span class="p">;</span></div><div class='line' id='LC590'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_third_jet_eta_ele</span><span class="p">;</span></div><div class='line' id='LC591'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_fourth_jet_eta_ele</span><span class="p">;</span></div><div class='line' id='LC592'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_ht1_ele</span><span class="p">;</span></div><div class='line' id='LC593'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_ht2_ele</span><span class="p">;</span></div><div class='line' id='LC594'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_ht3_ele</span><span class="p">;</span></div><div class='line' id='LC595'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_ht4_ele</span><span class="p">;</span></div><div class='line' id='LC596'><br/></div><div class='line' id='LC597'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_excmult_jets_muo</span><span class="p">;</span></div><div class='line' id='LC598'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_incmult_jets_muo</span><span class="p">;</span></div><div class='line' id='LC599'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_leading_jet_pt_muo</span><span class="p">;</span></div><div class='line' id='LC600'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_second_jet_pt_muo</span><span class="p">;</span></div><div class='line' id='LC601'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_third_jet_pt_muo</span><span class="p">;</span></div><div class='line' id='LC602'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_fourth_jet_pt_muo</span><span class="p">;</span></div><div class='line' id='LC603'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_leading_jet_eta_muo</span><span class="p">;</span></div><div class='line' id='LC604'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_second_jet_eta_muo</span><span class="p">;</span></div><div class='line' id='LC605'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_third_jet_eta_muo</span><span class="p">;</span></div><div class='line' id='LC606'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_fourth_jet_eta_muo</span><span class="p">;</span></div><div class='line' id='LC607'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_ht1_muo</span><span class="p">;</span></div><div class='line' id='LC608'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_ht2_muo</span><span class="p">;</span></div><div class='line' id='LC609'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_ht3_muo</span><span class="p">;</span></div><div class='line' id='LC610'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_ht4_muo</span><span class="p">;</span></div><div class='line' id='LC611'><br/></div><div class='line' id='LC612'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_excmult_jets_tot</span><span class="p">;</span></div><div class='line' id='LC613'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_incmult_jets_tot</span><span class="p">;</span></div><div class='line' id='LC614'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_leading_jet_pt_tot</span><span class="p">;</span></div><div class='line' id='LC615'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_second_jet_pt_tot</span><span class="p">;</span></div><div class='line' id='LC616'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_third_jet_pt_tot</span><span class="p">;</span></div><div class='line' id='LC617'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_fourth_jet_pt_tot</span><span class="p">;</span></div><div class='line' id='LC618'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_leading_jet_eta_tot</span><span class="p">;</span></div><div class='line' id='LC619'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_second_jet_eta_tot</span><span class="p">;</span></div><div class='line' id='LC620'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_third_jet_eta_tot</span><span class="p">;</span></div><div class='line' id='LC621'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_fourth_jet_eta_tot</span><span class="p">;</span></div><div class='line' id='LC622'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_ht1_tot</span><span class="p">;</span></div><div class='line' id='LC623'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_ht2_tot</span><span class="p">;</span></div><div class='line' id='LC624'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_ht3_tot</span><span class="p">;</span></div><div class='line' id='LC625'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="n">AIDA</span><span class="o">::</span><span class="n">IHistogram1D</span><span class="o">*</span> <span class="n">_h_ht4_tot</span><span class="p">;</span></div><div class='line' id='LC626'><br/></div><div class='line' id='LC627'>&nbsp;&nbsp;<span class="p">};</span></div><div class='line' id='LC628'><br/></div><div class='line' id='LC629'>&nbsp;&nbsp;<span class="c1">// This global object acts as a hook for the plugin system</span></div><div class='line' id='LC630'>&nbsp;&nbsp;<span class="n">AnalysisBuilder</span><span class="o">&lt;</span><span class="n">CMS_SMP_12_017</span><span class="o">&gt;</span> <span class="n">plugin_CMS_SMP_12_017</span><span class="p">;</span></div><div class='line' id='LC631'><br/></div><div class='line' id='LC632'><span class="p">}</span></div></pre></div></td>
         </tr>
       </table>
  </div>

  </div>
</div>

<a href="#jump-to-line" rel="facebox[.linejump]" data-hotkey="l" class="js-jump-to-line" style="display:none">Jump to Line</a>
<div id="jump-to-line" style="display:none">
  <form accept-charset="UTF-8" class="js-jump-to-line-form">
    <input class="linejump-input js-jump-to-line-field" type="text" placeholder="Jump to line&hellip;" autofocus>
    <button type="submit" class="button">Go</button>
  </form>
</div>

        </div>

      </div><!-- /.repo-container -->
      <div class="modal-backdrop"></div>
    </div><!-- /.container -->
  </div><!-- /.site -->


    </div><!-- /.wrapper -->

      <div class="container">
  <div class="site-footer">
    <ul class="site-footer-links right">
      <li><a href="https://status.github.com/">Status</a></li>
      <li><a href="http://developer.github.com">API</a></li>
      <li><a href="http://training.github.com">Training</a></li>
      <li><a href="http://shop.github.com">Shop</a></li>
      <li><a href="/blog">Blog</a></li>
      <li><a href="/about">About</a></li>

    </ul>

    <a href="/">
      <span class="mega-octicon octicon-mark-github" title="GitHub"></span>
    </a>

    <ul class="site-footer-links">
      <li>&copy; 2014 <span title="0.08742s from github-fe118-cp1-prd.iad.github.net">GitHub</span>, Inc.</li>
        <li><a href="/site/terms">Terms</a></li>
        <li><a href="/site/privacy">Privacy</a></li>
        <li><a href="/security">Security</a></li>
        <li><a href="/contact">Contact</a></li>
    </ul>
  </div><!-- /.site-footer -->
</div><!-- /.container -->


    <div class="fullscreen-overlay js-fullscreen-overlay" id="fullscreen_overlay">
  <div class="fullscreen-container js-fullscreen-container">
    <div class="textarea-wrap">
      <textarea name="fullscreen-contents" id="fullscreen-contents" class="fullscreen-contents js-fullscreen-contents" placeholder="" data-suggester="fullscreen_suggester"></textarea>
    </div>
  </div>
  <div class="fullscreen-sidebar">
    <a href="#" class="exit-fullscreen js-exit-fullscreen tooltipped tooltipped-w" aria-label="Exit Zen Mode">
      <span class="mega-octicon octicon-screen-normal"></span>
    </a>
    <a href="#" class="theme-switcher js-theme-switcher tooltipped tooltipped-w"
      aria-label="Switch themes">
      <span class="octicon octicon-color-mode"></span>
    </a>
  </div>
</div>



    <div id="ajax-error-message" class="flash flash-error">
      <span class="octicon octicon-alert"></span>
      <a href="#" class="octicon octicon-x close js-ajax-error-dismiss" aria-label="Dismiss error"></a>
      Something went wrong with that request. Please try again.
    </div>


      <script crossorigin="anonymous" src="https://assets-cdn.github.com/assets/frameworks-66f40787248d40d8457b357437c30ec844b57b28.js" type="text/javascript"></script>
      <script async="async" crossorigin="anonymous" src="https://assets-cdn.github.com/assets/github-266cb649e036760b856a495a9c114c49a592ae99.js" type="text/javascript"></script>
      
      
        <script async src="https://www.google-analytics.com/analytics.js"></script>
  </body>
</html>

