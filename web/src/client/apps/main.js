define(['knockout', 'crossroads', 'hasher'], function(ko, crossroads, hasher) {
    require('../components/index');

    var go_to_demo_default = function() {
        // HACK: This mantain the autoscroll functionality when reload.
        setTimeout(function() {
            if ($('[href="#demo"]')[0]) {
                // TODO: Detect page reload inside the component.
                $('[href="#demo"]')[0].click();
            }
        }, 500);
    };

    var go_to_demo = go_to_demo_default;

    var set_main_tag = function(tag) {
        // Jump to the demo area.
        go_to_demo();
        $('.main-zone').empty();
        $('.main-zone').append(tag);
        ko.applyBindings({}, tag);
    };

    var new_job_request = function() {
        go_to_demo = function() {};
        hasher.setHash('demo');
        go_to_demo = () => {
            go_to_demo = go_to_demo_default;
        };
    };

    // setup crossroads
    crossroads.addRoute('/demo', function() {
        var $jobRequest = $('<home></home>');
        set_main_tag($jobRequest[0]);
    });
    crossroads.addRoute('/help', function() {
        var $jobResult = $('<help></help>');
        set_main_tag($jobResult[0]);
    });
    crossroads.addRoute(':rest*:', new_job_request, -Infinity);

    // crossroads.routed.add(console.log, console); //log all routes

    // setup hasher
    function parseHash(newHash) {
        crossroads.parse(newHash);
    }
    hasher.initialized.add(parseHash); // parse initial hash
    hasher.changed.add(parseHash); // parse hash changes
    hasher.init(); // start listening for history change*/
});
