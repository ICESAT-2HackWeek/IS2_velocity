from IS2_velocity.readers import atl06_to_dict


def test_atl06_reader():
    data_dir = './data/'
    fn = 'processed_ATL06_20190822153035_08480411_003_01.h5'
    atl06_dict = atl06_to_dict(data_dir+fn,'/gt2l', index=None, epsg=3031)

    assert 'h_li' in atl06_dict
    assert 'x_atc' in atl06_dict
    assert 'x' in atl06_dict
