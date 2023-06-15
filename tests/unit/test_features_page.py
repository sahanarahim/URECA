NUM_CARDS, NUM_PICTURES = 6, 3

def test_features_page_should_be_accessible(client):
    response = client.get('/features')
    assert response.status_code == 200

def test_features_page_sections_should_render(client):
    response = client.get('/features')
    assert b'Features of our Database' in response.data
    assert b'Composition of PlantConnectome' in response.data
    assert b'Searching PlantConnectome' in response.data
    assert b'Search Results' in response.data

def test_piechart_canvas_should_render(client):
    response = client.get('/features')
    assert b"<canvas id = \"piechart\" width = '550' height = '200'> </canvas>" in response.data

def test_card_elements_should_render(client):
    response = client.get('/features')
    assert response.data.count(b"<div class = 'card'>") == NUM_CARDS
    assert response.data.count(b"<div class = 'card-divider'>") == NUM_CARDS
    assert response.data.count(b"<div class = 'card-section'>") == NUM_CARDS

def test_pictures_should_render(client):
    response = client.get('/features')
    assert response.data.count(b'<img') == NUM_PICTURES
