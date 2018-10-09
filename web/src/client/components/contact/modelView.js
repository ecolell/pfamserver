define('contact', [], function() {
    require('./style.css');

    function ContactModel(params) {
        this.init = function() {
            $('#contactForm input,#contactForm textarea').jqBootstrapValidation({
                preventSubmit: true,
                submitSuccess: function(not_used_variable, event) { // eslint-disable-line no-unused-vars
                    event.preventDefault(); // prevent default submit behaviour
                    // get values from FORM
                    var name = $('input#name').val();
                    var email = $('input#email').val();
                    var phone = $('input#phone').val();
                    var message = $('textarea#message').val();
                    var firstName = name; // For Success/Failure Message
                    // Check for white space in name for Success/Fail message
                    if (firstName.indexOf(' ') >= 0) {
                        firstName = name.split(' ').slice(0, -1).join(' ');
                    }
                    $('#contactForm button').html('<i class="fa icon-loading fa-spin"></i>  Sending Message');
                    $('input#name').prop('disabled', true);
                    $('input#email').prop('disabled', true);
                    $('input#phone').prop('disabled', true);
                    $('textarea#message').prop('disabled', true);
                    $('#contactForm button').prop('disabled', true);
                    $.ajax({
                        type: 'POST',
                        dataType: 'json',
                        contentType: 'application/json',
                        data: JSON.stringify({
                            name: name,
                            phone: phone,
                            email: email,
                            message: message,
                        }),
                        url: '/api/v0/notifications',
                        cache: false,
                        success: function() {
                            // Success message
                            $('#success').html("<div class='alert alert-success'>");
                            $('#success > .alert-success').html("<button type='button' class='close' data-dismiss='alert' aria-hidden='true'>&times;")
                                .append('</button>');
                            $('#success > .alert-success')
                                .append('<strong>Your message has been sent. </strong>');
                            $('#success > .alert-success')
                                .append('</div>');

                            // clear all fields
                            $('#contactForm').trigger('reset');
                            $('#contactForm button').html('Send Message');
                            $('input#name').prop('disabled', false);
                            $('input#email').prop('disabled', false);
                            $('input#phone').prop('disabled', false);
                            $('textarea#message').prop('disabled', false);
                            $('#contactForm button').prop('disabled', false);
                        },
                        error: function() {
                            // Fail message
                            $('#success').html("<div class='alert alert-danger'>");
                            $('#success > .alert-danger').html("<button type='button' class='close' data-dismiss='alert' aria-hidden='true'>&times;")
                                .append('</button>');
                            $('#success > .alert-danger').append($('<strong>').text('Sorry ' + firstName + ', it seems that my mail server is not responding. Please try again later!'));
                            $('#success > .alert-danger').append('</div>');
                            // clear all fields
                            $('#contactForm').trigger('reset');
                            $('#contactForm button').html('Send Message');
                            $('input#name').prop('disabled', false);
                            $('input#email').prop('disabled', false);
                            $('input#phone').prop('disabled', false);
                            $('textarea#message').prop('disabled', false);
                            $('#contactForm button').prop('disabled', false);
                        },
                    });
                },
                filter: function() {
                    return $(this).is(':visible');
                },
            });
        };
        this.init(params);
    }

    return ContactModel;
});
