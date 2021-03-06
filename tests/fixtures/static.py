import pytest
import re


@pytest.fixture(scope="function")
def main_app_template():
    return """<!DOCTYPE html>
<html lang="en">

<head>
    <link rel="shortcut icon" href="/static/favicon.ico">
    <!-- Google Tag Manager -->
    <script>(function(w,d,s,l,i){w[l]=w[l]||[];w[l].push({'gtm.start':
    new Date().getTime(),event:'gtm.js'});var f=d.getElementsByTagName(s)[0],
    j=d.createElement(s),dl=l!='dataLayer'?'&l='+l:'';j.async=true;j.src=
    'https://www.googletagmanager.com/gtm.js?id='+i+dl;f.parentNode.insertBefore(j,f);
    })(window,document,'script','dataLayer','GTM-KB5R6HC');</script>
    <!-- End Google Tag Manager -->

    <meta name="google-site-verification" content="GI8mhuV576ACMG70ypzNxitXHQwlyNQj8Twzu0v5mlU" />
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <meta name="description" content="">
    <meta name="author" content="">

    <title>pfamserver</title>

    <!-- HTML5 Shim and Respond.js IE8 support of HTML5 elements and media queries -->
    <!-- WARNING: Respond.js doesn't work if you view the page via file:// -->
    <!--[if lt IE 9]>


    <![endif]-->
    <link rel="stylesheet" href="http://localhost:2992/assets/main.2cd6ed.css">

</head>

<body id="page-top" class="index">
<!-- Google Tag Manager (noscript) -->
<!--<noscript><iframe src="https://www.googletagmanager.com/ns.html?id=GTM-KB5R6HC"
height="0" width="0" style="display:none;visibility:hidden"></iframe></noscript>-->
<!-- End Google Tag Manager (noscript) -->
    <script src="http://localhost:2992/assets/main.b65365e8.js"></script>

</body>
</html>""".split('\n')
