def test_entity_page_button_catalogue_should_render(client):
    response = client.get('/catalogue')
    assert b'<div class=\'button-group hollow\' style="padding-top: 10px; padding-bottom: 10px;">' in response.data

def test_entity_page_table_and_pagination_should_render(client):
    response = client.get('/catalogue')
    assert b'<table>' in response.data
    assert b'<ul class="pagination text-center">' in response.data