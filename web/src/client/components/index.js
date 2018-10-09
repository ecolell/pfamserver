define(['require', 'knockout', 'mini-toastr', 'lodash/forEach'], function(require, ko, miniToastr, forEach) {
    // services
    miniToastr.init();

    // components

    function register(name) {
        ko.components.register(name, {
            viewModel: require('./' + name + '/modelView.js'),
            template: require('./' + name + '/view.html'),
        });
    }

    var components = [
        'main',
        'home',
        'examples',
        'team',
        'contact',
    ];

    forEach(components, register);
});
