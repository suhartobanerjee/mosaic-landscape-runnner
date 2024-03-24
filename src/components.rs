use bevy::window::PrimaryWindow;
use bevy::prelude::*;

pub const PLAYER_SPEED: f32 = 500.0;
pub const BASE_ROAD: (f32, f32, f32)= (0.0, 70.0, 0.0);
pub const BASE_PLAYER: (f32, f32, f32)= (80.0, BASE_ROAD.1 + 60.0, BASE_ROAD.2);
pub const PLAYER_SIZE: f32 = 70.0;


#[derive(Component)]
pub struct Player {}

pub fn spawn_player(mut commands: Commands,
                    window_query: Query<&Window, With<PrimaryWindow>>,
                    asset_server: Res<AssetServer>) {
    
    //let window = window_query.get_single().unwrap();

    commands.spawn((
        SpriteBundle{
            transform: Transform::from_xyz(BASE_PLAYER.0, BASE_PLAYER.1, BASE_PLAYER.2),
            texture: asset_server.load("sprites/character_maleAdventurer_run0.png"),
            ..Default::default()
        },
        Player {},
        )
    );
}

pub fn spawn_camera(mut commands: Commands,
                    window_query: Query<&Window, With<PrimaryWindow>>) {
    
    let window = window_query.get_single().unwrap();

    commands.spawn(
        Camera2dBundle {
            transform: Transform::from_xyz(window.width() / 2.0, window.height() / 2.0, BASE_ROAD.2),
            ..Default::default()
        }
    );
}


pub fn player_movement(mut commands: Commands,
                       keyboard_input: Res<Input<KeyCode>>,
                       mut player_query: Query<&mut Transform, With<Player>>,
                       window_query: Query<&Window, With<PrimaryWindow>>,
                       time: Res<Time>) {

    if let Ok(mut transform) = player_query.get_single_mut() {
        let mut direction = Vec3::ZERO;
        let window = window_query.get_single().unwrap();

        if keyboard_input.pressed(KeyCode::Right) {
            direction += Vec3::new(1.0, 0.0, 0.0)
        }
        if keyboard_input.pressed(KeyCode::Left) {
            direction += Vec3::new(-1.0, 0.0, 0.0)
        }
        if keyboard_input.pressed(KeyCode::Up) {
            direction += Vec3::new(0.0, 1.0, 0.0)
        }
        if keyboard_input.pressed(KeyCode::Down) {
            direction += Vec3::new(0.0, -1.0, 0.0)
        }

        if direction.length() > 0.0 {
            direction = direction.normalize();
        }

        transform.translation += direction * PLAYER_SPEED * time.delta_seconds();
        commands.spawn(
            Camera2dBundle {
                transform: Transform::from_xyz(transform.translation.x, transform.translation.y, transform.translation.z),
                ..Default::default()
            }
            );
    }
}


pub fn confine_player_movement(mut player_query: Query<&mut Transform, With<Player>>,
                               window_query: Query<&Window, With<PrimaryWindow>>) {

    if let Ok(mut player_transform) = player_query.get_single_mut() {
        
        let window = window_query.get_single().unwrap();
        let x_min = BASE_PLAYER.0;
        let x_max = window.width() - PLAYER_SIZE / 2.0;
        let y_min = BASE_PLAYER.1;
        let y_max = y_min + BASE_ROAD.1;

        let mut translation = player_transform.translation;

        //bound the movement
        if translation.x < x_min {
           translation.x = x_min; 
        }
        if translation.x > x_max {
           translation.x = x_max; 
        }
        if translation.y > y_max {
           translation.y = y_max; 
        }
        if translation.y < y_min {
           translation.y = y_min; 
        }

        player_transform.translation = translation;
    }
    
}


#[derive(Component)]
pub struct Landscape {
    sv_state: Vec<String>
}

pub fn spawn_landscape(mut commands: Commands,
                       //window_query: Query<&Window, With<PrimaryWindow>>,
                       asset_server: Res<AssetServer>) {
    
    //let window = window_query.get_single().unwrap();

    let mut texture_atlas = TextureAtlasBuilder::default();
    //texture_atlas.add_texture(asset_server.load("sprites/roadTexture_13.png"), ..default());

    let offset = 120.0;
    let rounds: u16 = 12;
    for round in 0..rounds{

        commands.spawn((
                SpriteBundle{
                    transform: Transform::from_xyz(BASE_ROAD.0 + offset * f32::from(round), BASE_ROAD.1, BASE_ROAD.2),
                    texture: asset_server.load("sprites/roadTexture_13.png"),
                    ..Default::default()
                },
                Landscape {
                    sv_state: vec!["ref".to_string(), "ref".to_string()]
                },
                )
                      );
    }
}


