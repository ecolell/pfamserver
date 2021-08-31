define('home', ['hasher', 'knockout', 'lodash/forEach'], function(hasher, ko, forEach) {
    require('bootstrap/dist/css/bootstrap.css');
    require('./style.font-awesome.css');
    require('url-loader!./style.fonts.css');
    require('./style.agency.css');
    require('./style.leloir.css');
    require('./style.dropzone.css');

    function HomeModel(params) {
        var self = this;

        this.request = params.request;
        this.public_key = params.public_key;

        this.demo_label = ko.computed(function() {
            return self.request ? 'Start your job' : 'Your findings';
        });

        this.demo_subheading = ko.computed(function() {
            return self.request ? 'Tell us more about your protein.' : 'Take a look at the results!';
        });

        this.openHelp = function() {
            hasher.setHash('help');
        };

        this.tell_me_more = function() {
            $('html, body').animate({ scrollTop: $('#services').offset().top }, 1000);
        };

        this.start_new_job = function() {
            if (self.request) {
                $('html, body').animate({ scrollTop: $('#demo').offset().top }, 1000);
            } else {
                hasher.setHash('/job_request/~');
            }
        };

        this.initializeAgency = function() {
            $('a.page-scroll').bind('click', function(event) {
                var $anchor = $(this);
                $('html, body').stop().animate({
                    scrollTop: ($($anchor.attr('href')).offset().top - 50),
                }, 1250, 'easeInOutExpo');
                event.preventDefault();
            });

            $('body').scrollspy({
                target: '.navbar-fixed-top',
                offset: 51,
            });

            // Closes the Responsive Menu on Menu Item Click
            $('.navbar-collapse ul li a').click(function() {
                $('.navbar-toggle:visible').click();
            });

            $('#mainNav').affix({
                offset: {
                    top: 100,
                },
            });
        };

        this.initializeContactUS = function() {
            $('a[data-toggle="tab"]').click(function(e) {
                e.preventDefault();
                $(this).tab('show');
            });

            /* When clicking on Full hide fail/success boxes */
            $('#name').focus(function() {
                $('#success').html('');
            });
        };

        this.init = function() {
            forEach($('img[data-src]'), function(img) {
                img.setAttribute('src', img.getAttribute('data-src'));
            });

            self.initializeAgency();
            self.initializeContactUS();
        };
        this.init(params);
    }

    return HomeModel;
});
