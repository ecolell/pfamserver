define('main', [], function() {
    var body = require('./view.html');
    $('body').append(body);

    // Plugin JavaScript
    require('jquery-easing');
    // Bootstrap Core JavaScript
    require('bootstrap/dist/js/bootstrap.js');
    // Contact Form JavaScript
    require('./vendor/jqBootstrapValidation.js');

    function MainModel() {

    }

    return MainModel;
});
