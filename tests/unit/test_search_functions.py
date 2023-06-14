def test_search_routes_should_only_accept_get_requests(client):
    pass
    '''
    Fill this out later
    '''
    for i in ['normal', 'exact', 'alias', 'substring', 'non_alpha']:
        req_post, req_put = client.post(), client.put()
        req_patch, req_del = client.patch(), client.delete()
        
        assert req_post.status_code == 405
        assert req_put.status_code == 405
        assert req_patch.status_code == 405
        assert req_del.status_code == 405